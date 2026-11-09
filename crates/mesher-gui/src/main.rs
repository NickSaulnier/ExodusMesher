//! `mesher-gui`: pick parts of an OBJ file, run the worker, export. It holds
//! no meshing logic: models are read and diagnosed with the worker's own code
//! for display (on a background thread), and meshing runs in `mesher-worker`
//! through `mesher-client`, exactly as the CLI does. Every export saves a
//! project file that `mesher run project.toml` replays.

// Camera and display maths only: nothing here reaches a mesh or a report, so
// the platform's maths library is fine (see `clippy.toml`).
#![allow(clippy::disallowed_methods)]
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod render;
mod theme;

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::time::Duration;

use eframe::egui::{self, Color32, RichText};
use mesher_core::project::Project;
use mesher_core::protocol::{MESH_EXODUS_FILE, WorkerMessage};
use mesher_core::{ExitStatus, Record, Severity};
use mesher_gui::camera::Camera;
use mesher_gui::scene::{Scene, severity_slot};
use mesher_gui::session::Session;
use render::{Uniforms, Vertex, Viewport};
use theme::{Icon, IconButton};

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("ExodusMesher")
            .with_inner_size([1440.0, 900.0])
            .with_min_inner_size([960.0, 600.0]),
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };
    let start = std::env::args_os().nth(1).map(PathBuf::from);
    eframe::run_native(
        "ExodusMesher",
        options,
        Box::new(move |cc| Ok(Box::new(App::new(cc, start)))),
    )
}

/// Part colours (sRGB); unselected parts are grey.
const PALETTE: [[u8; 3]; 8] = [
    [76, 140, 220],
    [230, 150, 60],
    [90, 180, 110],
    [200, 90, 160],
    [150, 120, 210],
    [70, 180, 190],
    [210, 190, 80],
    [190, 110, 90],
];
const UNSELECTED: [u8; 3] = [150, 150, 150];
const DEFECT: [u8; 3] = [235, 60, 50];
const HIGHLIGHT: [u8; 3] = [255, 220, 40];

/// A colour for the viewport: sRGB components in [0, 1] (egui shows the
/// viewport texture's values as they are, so no linearisation).
fn rgba(c: [u8; 3], alpha: f32) -> [f32; 4] {
    let l = |v: u8| f32::from(v) / 255.0;
    [l(c[0]), l(c[1]), l(c[2]), alpha]
}

fn part_colour(i: usize) -> [u8; 3] {
    PALETTE[i % PALETTE.len()]
}

fn severity_colour(s: Severity) -> Color32 {
    match s {
        Severity::Info => Color32::from_rgb(120, 160, 220),
        Severity::Warning => Color32::from_rgb(220, 170, 50),
        Severity::Error => Color32::from_rgb(230, 100, 60),
        Severity::Fatal => Color32::from_rgb(230, 50, 50),
    }
}

/// A background job's news.
enum RunEvent {
    Message(WorkerMessage),
    Done { exit_code: i32, summary: String },
}

struct Run {
    rx: Receiver<RunEvent>,
    cancel: Arc<AtomicBool>,
    stage: String,
    fraction: f32,
}

struct App {
    session: Option<Session>,
    loading: Option<(PathBuf, Receiver<Scene>, Option<Project>)>,
    camera: Camera,
    viewport: Option<Viewport>,
    dirty: bool,
    wireframe: bool,
    defects: bool,
    clip: bool,
    clip_at: f32,
    record: Option<usize>,
    /// Which severities the records list shows, by `severity_slot`.
    record_filter: [bool; 4],
    run: Option<Run>,
    log: Vec<(Severity, String)>,
    status: String,
    worker: PathBuf,
    timeout_seconds: f64,
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>, start: Option<PathBuf>) -> App {
        theme::apply(&cc.egui_ctx);
        let viewport = cc.wgpu_render_state.as_ref().map(Viewport::new);
        let mut app = App {
            session: None,
            loading: None,
            camera: Camera::fit([[0.0; 3], [1.0; 3]]),
            viewport,
            dirty: true,
            wireframe: true,
            defects: true,
            clip: false,
            clip_at: 0.5,
            record: None,
            record_filter: [true; 4],
            run: None,
            log: Vec::new(),
            status: "Open an OBJ file or a project to begin.".to_owned(),
            worker: mesher_client::default_worker_path(),
            timeout_seconds: 3600.0,
        };
        if let Some(path) = start {
            app.open(&path, &cc.egui_ctx);
        }
        app
    }

    /// Open an OBJ file, or a project file (its input, selection and
    /// settings), reading it on a background thread.
    fn open(&mut self, path: &Path, ctx: &egui::Context) {
        let is_project = path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("toml"));
        let (input, project) = if is_project {
            match Project::read(path) {
                Ok(p) => (p.input.clone(), Some(p)),
                Err(e) => {
                    self.status = format!("Cannot open project: {e}");
                    return;
                }
            }
        } else {
            (path.to_path_buf(), None)
        };
        let (tx, rx) = mpsc::channel();
        let ctx = ctx.clone();
        let read = input.clone();
        std::thread::spawn(move || {
            let _ = tx.send(Scene::load(&read));
            ctx.request_repaint();
        });
        self.status = format!("Reading {} ...", short_path(&input));
        self.loading = Some((input, rx, project));
    }

    fn poll_loading(&mut self) {
        let Some((path, rx, _)) = &self.loading else {
            return;
        };
        match rx.try_recv() {
            Ok(scene) => {
                let project = self.loading.take().and_then(|(_, _, p)| p);
                self.camera = Camera::fit(scene.bounds);
                self.clip_at = 0.5;
                self.record = None;
                self.status = format!(
                    "{}: {} parts, {} triangles, {} records",
                    short_path(&scene.path),
                    scene.parts.len(),
                    scene.model.triangles.len(),
                    scene.records.len()
                );
                let mut session = Session::new(scene);
                if let Some(p) = &project {
                    session.apply(p);
                }
                self.session = Some(session);
                self.dirty = true;
            }
            Err(TryRecvError::Empty) => {}
            Err(TryRecvError::Disconnected) => {
                self.status = format!("Reading {} failed", short_path(path));
                self.loading = None;
            }
        }
    }

    /// Start the worker on the session's project: create the output
    /// directory, save the project beside it, and stream progress.
    fn start_run(&mut self, ctx: &egui::Context) {
        let Some(session) = &self.session else {
            return;
        };
        if !session.any_selected() {
            self.status = "Select at least one part.".to_owned();
            return;
        }
        let project = session.project();
        if let Err(e) = std::fs::create_dir_all(&session.output) {
            self.status = format!("Cannot create {}: {e}", short_path(&session.output));
            return;
        }
        let project_path = session.project_path();
        if let Err(e) = project.write(&project_path) {
            self.status = format!("Cannot save the project: {e}");
            return;
        }
        self.log.clear();
        self.log.push((
            Severity::Info,
            format!(
                "Saved {} (replay: mesher run \"{}\")",
                project_path.display(),
                project_path.display()
            ),
        ));
        let (tx, rx) = mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let stop = Arc::clone(&cancel);
        let worker = self.worker.clone();
        let limits = mesher_client::Limits {
            timeout: Duration::try_from_secs_f64(self.timeout_seconds)
                .unwrap_or(Duration::from_secs(3600)),
            memory: Some(8192 << 20),
        };
        let job = project.job();
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let mut forward = |m: &WorkerMessage| {
                let _ = tx.send(RunEvent::Message(m.clone()));
                ctx.request_repaint();
            };
            let (report, _) = mesher_client::run_job(&worker, &job, limits, &mut forward, &stop);
            let s = &report.summary;
            let summary = format!(
                "exit {}: {} info, {} warnings, {} errors, {} fatal",
                s.exit_code, s.info, s.warnings, s.errors, s.fatal
            );
            let _ = tx.send(RunEvent::Done {
                exit_code: s.exit_code,
                summary,
            });
            ctx.request_repaint();
        });
        self.status = "Running ...".to_owned();
        self.run = Some(Run {
            rx,
            cancel,
            stage: "starting".to_owned(),
            fraction: 0.0,
        });
    }

    fn poll_run(&mut self) {
        let Some(run) = &mut self.run else {
            return;
        };
        loop {
            match run.rx.try_recv() {
                Ok(RunEvent::Message(m)) => match m {
                    WorkerMessage::Progress { stage, fraction } => {
                        run.stage = format!("{stage:?}");
                        run.fraction = fraction as f32;
                    }
                    WorkerMessage::Record { record } => {
                        self.log.push((record.severity, record_line(&record)));
                    }
                    _ => {}
                },
                Ok(RunEvent::Done { exit_code, summary }) => {
                    let outcome = match ExitStatus::from_code(exit_code) {
                        Some(ExitStatus::Clean) => "clean",
                        Some(ExitStatus::Repaired) => "succeeded after repairs",
                        Some(ExitStatus::Fallback) => "succeeded via fallback",
                        Some(ExitStatus::Failed) => "failed; see the report",
                        Some(ExitStatus::InternalBug) => "internal error",
                        None => "unknown",
                    };
                    let out = self
                        .session
                        .as_ref()
                        .map(|s| short_path(&s.output.join(MESH_EXODUS_FILE)))
                        .unwrap_or_default();
                    self.status = format!("Run {outcome} ({summary}); Exodus: {out}");
                    self.log.push((Severity::Info, self.status.clone()));
                    self.run = None;
                    return;
                }
                Err(TryRecvError::Empty) => return,
                Err(TryRecvError::Disconnected) => {
                    self.status = "The run ended without a result".to_owned();
                    self.run = None;
                    return;
                }
            }
        }
    }

    /// Rebuild the viewport geometry from the session.
    fn rebuild(&mut self, frame: &eframe::Frame) {
        let (Some(rs), Some(viewport)) = (frame.wgpu_render_state(), self.viewport.as_mut()) else {
            return;
        };
        let (tris, lines) = match &self.session {
            Some(s) => geometry(s, self.wireframe, self.defects, self.record),
            None => (Vec::new(), Vec::new()),
        };
        viewport.set_geometry(rs, &tris, &lines);
        self.dirty = false;
    }

    fn parts_panel(&mut self, ui: &mut egui::Ui) {
        let Some(session) = &mut self.session else {
            ui.heading("Parts");
            ui.label(theme::dim("No model loaded."));
            return;
        };
        let chosen = session.selected.iter().filter(|&&s| s).count();
        ui.horizontal(|ui| {
            ui.heading("Parts");
            ui.label(theme::dim(format!(
                "{chosen} of {}",
                session.selected.len()
            )));
        });
        ui.horizontal(|ui| {
            if ui.add(IconButton::new(Icon::Check, "All")).clicked() {
                session.selected.iter_mut().for_each(|s| *s = true);
                self.dirty = true;
            }
            if ui.add(IconButton::new(Icon::Cross, "None")).clicked() {
                session.selected.iter_mut().for_each(|s| *s = false);
                self.dirty = true;
            }
        });
        ui.add_space(4.0);
        egui::ScrollArea::vertical().show(ui, |ui| {
            for i in 0..session.scene.parts.len() {
                if part_card(ui, session, i) {
                    session.toggle(i);
                    self.dirty = true;
                }
            }
        });
    }

    fn run_panel(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        ui.heading("Run");
        let running = self.run.is_some();
        let mut go = false;
        if let Some(session) = &mut self.session {
            let s = &mut session.settings;
            egui::Grid::new("settings")
                .num_columns(3)
                .spacing([14.0, 10.0])
                .show(ui, |ui| {
                    ui.label("Element size").on_hover_text(
                        "Target element size, relative to each part's bounding-box diagonal",
                    );
                    ui.add(percent(&mut s.element_size, 0.01, 100.0));
                    ui.label(theme::dim("of diagonal").small());
                    ui.end_row();
                    ui.label("Max deviation").on_hover_text(
                        "Largest acceptable fallback deviation, relative to the diagonal",
                    );
                    ui.add(percent(&mut s.max_deviation, 0.0001, 100.0));
                    ui.label(theme::dim("of diagonal").small());
                    ui.end_row();
                    ui.label("Strict");
                    theme::switch(ui, &mut s.strict);
                    ui.label(theme::dim("conformal only, no fallback").small());
                    ui.end_row();
                    ui.label("No repair");
                    theme::switch(ui, &mut s.no_repair);
                    ui.label(theme::dim("diagnose and report only").small());
                    ui.end_row();
                    ui.label("Uniform size");
                    theme::switch(ui, &mut s.uniform_size);
                    ui.label(theme::dim("no grading for curvature or thin walls").small());
                    ui.end_row();
                    ui.label("No remesh");
                    theme::switch(ui, &mut s.no_remesh);
                    ui.label(theme::dim("no remesh-and-retry before the fallback").small());
                    ui.end_row();
                });
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.label("Output");
                if ui.add(IconButton::new(Icon::Folder, "Choose...")).clicked() {
                    if let Some(dir) = rfd::FileDialog::new()
                        .set_directory(&session.output)
                        .pick_folder()
                    {
                        session.output = dir;
                    }
                }
            });
            egui::Frame::new()
                .fill(theme::SUNKEN)
                .corner_radius(6)
                .inner_margin(egui::Margin::symmetric(10, 6))
                .show(ui, |ui| {
                    ui.set_min_width(ui.available_width());
                    ui.add(
                        egui::Label::new(theme::dim(short_path(&session.output)).small())
                            .truncate(),
                    )
                    .on_hover_text(session.output.display().to_string());
                });
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.add_enabled_ui(!running && session.any_selected(), |ui| {
                    if ui
                        .add(IconButton::primary(Icon::Play, "Run and export"))
                        .on_hover_text("Mesh the selected parts and write Exodus, .msh, .vtu, the report and project.toml")
                        .clicked()
                    {
                        go = true;
                    }
                });
            });
        } else {
            ui.label("No model loaded.");
        }
        if go {
            self.start_run(ctx);
        }
        if let Some(run) = &self.run {
            ui.add(egui::ProgressBar::new(run.fraction).text(format!(
                "{} {:.0}%",
                run.stage,
                run.fraction * 100.0
            )));
            if ui.add(IconButton::new(Icon::Stop, "Cancel")).clicked() {
                run.cancel.store(true, Ordering::Relaxed);
            }
        }
        ui.separator();
        ui.collapsing("Run log", |ui| {
            egui::ScrollArea::vertical()
                .id_salt("log")
                .max_height(160.0)
                .stick_to_bottom(true)
                .show(ui, |ui| {
                    for (sev, line) in &self.log {
                        ui.label(RichText::new(line).small().color(severity_colour(*sev)));
                    }
                });
        });
        ui.separator();
        self.records_panel(ui);
    }

    fn records_panel(&mut self, ui: &mut egui::Ui) {
        let Some(session) = &self.session else {
            return;
        };
        let records = &session.scene.records;
        ui.horizontal(|ui| {
            ui.heading("Records");
            ui.label(theme::dim(records.len().to_string()));
        });
        // Severity filters, with counts.
        let mut counts = [0usize; 4];
        for r in records {
            counts[severity_slot(r.severity)] += 1;
        }
        ui.horizontal_wrapped(|ui| {
            for sev in [
                Severity::Info,
                Severity::Warning,
                Severity::Error,
                Severity::Fatal,
            ] {
                let slot = severity_slot(sev);
                if counts[slot] == 0 {
                    continue;
                }
                let on = self.record_filter[slot];
                let text = format!("{} {}", short(sev), counts[slot]);
                if theme::chip(ui, &text, severity_colour(sev), on)
                    .on_hover_text("Show or hide this severity")
                    .clicked()
                {
                    self.record_filter[slot] = !on;
                }
            }
        });
        if records.is_empty() {
            ui.label(theme::dim("No defects found.").small());
            return;
        }
        let mut clicked = None;
        let list_height = if self.record.is_some() { 220.0 } else { 420.0 };
        egui::ScrollArea::vertical()
            .id_salt("records")
            .max_height(list_height)
            .auto_shrink([false, true])
            .show(ui, |ui| {
                for (i, r) in records.iter().enumerate() {
                    if !self.record_filter[severity_slot(r.severity)] {
                        continue;
                    }
                    if record_row(ui, r, i, self.record == Some(i)) {
                        clicked = Some(i);
                    }
                }
            });
        if let Some(i) = clicked {
            self.record = if self.record == Some(i) {
                None
            } else {
                Some(i)
            };
            self.dirty = true;
        }
        if let Some(r) = self.record.and_then(|i| records.get(i)) {
            ui.add_space(6.0);
            record_details(ui, r);
        }
    }

    fn viewport_panel(&mut self, ui: &mut egui::Ui, frame: &eframe::Frame) {
        let rect = ui.available_rect_before_wrap();
        let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());
        let aspect = rect.width() / rect.height().max(1.0);
        if response.dragged_by(egui::PointerButton::Primary) {
            let d = response.drag_delta();
            self.camera.orbit(d.x, d.y);
        }
        if response.dragged_by(egui::PointerButton::Secondary)
            || response.dragged_by(egui::PointerButton::Middle)
        {
            let d = response.drag_delta();
            self.camera.pan(d.x, d.y, rect.height());
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                self.camera.zoom((-scroll * 0.002).exp());
            }
        }
        if response.clicked() {
            if let (Some(pos), Some(session)) = (response.interact_pointer_pos(), &mut self.session)
            {
                let x = (pos.x - rect.left()) / rect.width() * 2.0 - 1.0;
                let y = 1.0 - (pos.y - rect.top()) / rect.height() * 2.0;
                let (origin, dir) = self.camera.ray(x, y, aspect);
                if let Some(part) = session.scene.pick(origin, dir) {
                    session.toggle(part);
                    self.dirty = true;
                }
            }
        }
        if response.double_clicked() {
            if let Some(s) = &self.session {
                self.camera = Camera::fit(s.scene.bounds);
            }
        }

        if self.dirty {
            self.rebuild(frame);
        }
        let (Some(rs), Some(viewport)) = (frame.wgpu_render_state(), self.viewport.as_mut()) else {
            ui.painter().text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "The viewport needs the wgpu renderer",
                egui::FontId::proportional(18.0),
                Color32::GRAY,
            );
            return;
        };
        let ppp = ui.ctx().pixels_per_point();
        let size = [
            (rect.width() * ppp).round() as u32,
            (rect.height() * ppp).round() as u32,
        ];
        let bounds = self
            .session
            .as_ref()
            .map_or([[0.0; 3], [1.0; 3]], |s| s.scene.bounds);
        // Clip along x at a fraction of the bounding box.
        let at = bounds[0][0] + (bounds[1][0] - bounds[0][0]) * f64::from(self.clip_at);
        let back = self.camera.back();
        let uniforms = Uniforms {
            view_proj: self.camera.view_proj(aspect),
            light: [back[0], back[1], back[2], 0.0],
            clip: [1.0, 0.0, 0.0, at as f32],
            flags: [if self.clip { 1.0 } else { 0.0 }, 0.0, 0.0, 0.0],
        };
        let id = viewport.render(rs, size, &uniforms);
        ui.painter().image(
            id,
            rect,
            egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
            Color32::WHITE,
        );
        ui.painter().text(
            rect.left_bottom() + egui::vec2(12.0, -10.0),
            egui::Align2::LEFT_BOTTOM,
            "drag: orbit   right-drag: pan   scroll: zoom   click: toggle part   double-click: fit",
            egui::FontId::proportional(14.0),
            Color32::from_gray(170),
        );
    }
}

/// One part as a card: checkbox, colour swatch, name, size and defect
/// badges. Selected parts get a filled card; the whole card highlights under
/// the pointer and toggles on click. Returns whether it was toggled.
fn part_card(ui: &mut egui::Ui, session: &mut Session, i: usize) -> bool {
    let part = &session.scene.parts[i];
    let selected = session.selected[i];
    let [r, g, b] = part_colour(i);
    let size = {
        let [lo, hi] = part.bounds;
        (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
    };
    let mut frame = egui::Frame::new()
        .inner_margin(egui::Margin::symmetric(10, 8))
        .corner_radius(8)
        .begin(ui);
    let mut toggled = false;
    {
        let ui = &mut frame.content_ui;
        ui.set_min_width(ui.available_width());
        ui.horizontal(|ui| {
            let mut on = selected;
            toggled |= ui.checkbox(&mut on, "").changed();
            let (swatch, _) = ui.allocate_exact_size(egui::vec2(14.0, 14.0), egui::Sense::hover());
            let swatch_colour = if selected {
                Color32::from_rgb(r, g, b)
            } else {
                Color32::from_gray(90)
            };
            ui.painter().rect_filled(swatch, 4.0, swatch_colour);
            let name = RichText::new(&part.name).strong();
            ui.label(if selected {
                name
            } else {
                name.color(theme::TEXT_DIM)
            });
        });
        ui.horizontal(|ui| {
            ui.add_space(ui.spacing().icon_width + ui.spacing().item_spacing.x);
            ui.label(
                theme::dim(format!(
                    "{} triangles · size {size:.3}",
                    part.triangles.len()
                ))
                .small(),
            );
            for sev in [Severity::Warning, Severity::Error, Severity::Fatal] {
                let n = part.badges[severity_slot(sev)];
                if n > 0 {
                    ui.label(
                        RichText::new(format!("{} {n}", short(sev)))
                            .small()
                            .color(severity_colour(sev)),
                    );
                }
            }
        });
    }
    let outer = frame.allocate_space(ui);
    let card = ui.interact(
        outer.rect,
        egui::Id::unique(("part card", i)),
        egui::Sense::click(),
    );
    frame.frame.fill = if card.hovered() {
        theme::CARD_HOVER
    } else if selected {
        theme::CARD
    } else {
        Color32::TRANSPARENT
    };
    frame.paint(ui);
    toggled || card.clicked()
}

/// One record as a selectable row: severity chip, code and part, with the
/// hint below. Returns whether it was clicked.
fn record_row(ui: &mut egui::Ui, r: &Record, i: usize, selected: bool) -> bool {
    let mut frame = egui::Frame::new()
        .inner_margin(egui::Margin::symmetric(8, 6))
        .corner_radius(6)
        .begin(ui);
    {
        let ui = &mut frame.content_ui;
        ui.set_min_width(ui.available_width());
        ui.spacing_mut().item_spacing.y = 2.0;
        ui.spacing_mut().interact_size.y = 22.0;
        ui.horizontal(|ui| {
            theme::chip(ui, short(r.severity), severity_colour(r.severity), true);
            ui.label(RichText::new(r.code.as_str()).strong().small());
            if let Some(part) = &r.part {
                ui.label(theme::dim(part.as_str()).small());
            }
        });
        if let Some(hint) = &r.hint {
            ui.add(egui::Label::new(theme::dim(hint.as_str()).small()).truncate());
        }
    }
    let outer = frame.allocate_space(ui);
    let row = ui.interact(
        outer.rect,
        egui::Id::unique(("record row", i)),
        egui::Sense::click(),
    );
    frame.frame.fill = if row.hovered() {
        theme::CARD_HOVER
    } else if selected {
        theme::CARD
    } else {
        Color32::TRANSPARENT
    };
    if selected {
        frame.frame.stroke = egui::Stroke::new(1.0, theme::ACCENT);
    }
    frame.paint(ui);
    row.clicked()
}

/// The chosen record in full: description, hint, measurements against their
/// thresholds, and where it is in the OBJ file.
fn record_details(ui: &mut egui::Ui, r: &Record) {
    egui::Frame::new()
        .fill(theme::SUNKEN)
        .corner_radius(8)
        .inner_margin(egui::Margin::same(12))
        .show(ui, |ui| {
            ui.set_min_width(ui.available_width());
            ui.spacing_mut().interact_size.y = 20.0;
            ui.horizontal(|ui| {
                theme::chip(ui, short(r.severity), severity_colour(r.severity), true);
                ui.label(RichText::new(r.code.as_str()).strong());
            });
            ui.label(r.code.description());
            if let Some(hint) = &r.hint {
                ui.label(theme::dim(hint.as_str()));
            }
            ui.add_space(4.0);
            egui::Grid::new("record details")
                .num_columns(2)
                .spacing([12.0, 4.0])
                .show(ui, |ui| {
                    if let Some(part) = &r.part {
                        ui.label(theme::dim("part").small());
                        ui.label(RichText::new(part).small());
                        ui.end_row();
                    }
                    for (k, v) in &r.measured {
                        ui.label(theme::dim(k.as_str()).small());
                        ui.label(RichText::new(value_text(v)).small().monospace());
                        ui.end_row();
                    }
                    for (k, v) in &r.threshold {
                        ui.label(theme::dim(format!("{k} (limit)")).small());
                        ui.label(RichText::new(value_text(v)).small().monospace());
                        ui.end_row();
                    }
                    let lines = &r.location.obj_lines;
                    if !lines.is_empty() {
                        let shown: Vec<String> =
                            lines.iter().take(12).map(u32::to_string).collect();
                        ui.label(theme::dim("OBJ lines").small());
                        ui.label(
                            RichText::new(format!(
                                "{}{}",
                                shown.join(", "),
                                if lines.len() > 12 { ", ..." } else { "" }
                            ))
                            .small()
                            .monospace(),
                        );
                        ui.end_row();
                    }
                });
        });
}

/// A measured or threshold value for display: strings without quotes,
/// floats to about four significant digits.
fn value_text(v: &serde_json::Value) -> String {
    if let Some(s) = v.as_str() {
        return s.to_owned();
    }
    match v.as_f64() {
        Some(x) if v.is_f64() => {
            if x == 0.0 || (1e-3..1e5).contains(&x.abs()) {
                let digits = (3 - x.abs().log10().floor() as i32).clamp(0, 6) as usize;
                format!("{x:.digits$}")
            } else {
                format!("{x:.3e}")
            }
        }
        _ => v.to_string(),
    }
}

/// A drag field for a fraction, shown and typed as a percentage.
fn percent(value: &mut f64, min_percent: f64, max_percent: f64) -> egui::DragValue<'_> {
    egui::DragValue::new(value)
        .speed(0.0002)
        .range(min_percent / 100.0..=max_percent / 100.0)
        .custom_formatter(|v, _| {
            let p = v * 100.0;
            if p >= 1.0 {
                format!("{p:.2} %")
            } else {
                format!("{p:.3} %")
            }
        })
        .custom_parser(|text| {
            text.trim()
                .trim_end_matches('%')
                .trim()
                .parse::<f64>()
                .ok()
                .map(|p| p / 100.0)
        })
}

/// A path shortened for display: the full path when short, else its last
/// two components after an ellipsis.
fn short_path(path: &Path) -> String {
    let full = path.display().to_string();
    if full.chars().count() <= 48 {
        return full;
    }
    let mut tail: Vec<String> = path
        .components()
        .rev()
        .take(2)
        .map(|c| c.as_os_str().to_string_lossy().into_owned())
        .collect();
    tail.reverse();
    let sep = std::path::MAIN_SEPARATOR.to_string();
    format!("...{sep}{}", tail.join(&sep))
}

fn short(s: Severity) -> &'static str {
    match s {
        Severity::Info => "info",
        Severity::Warning => "warn",
        Severity::Error => "error",
        Severity::Fatal => "fatal",
    }
}

fn record_line(r: &Record) -> String {
    let part = r
        .part
        .as_deref()
        .map(|p| format!(" [{p}]"))
        .unwrap_or_default();
    let hint = r
        .hint
        .as_deref()
        .map(|h| format!(": {h}"))
        .unwrap_or_default();
    format!("{} {}{part}{hint}", short(r.severity), r.code.as_str())
}

/// Viewport geometry: flat-shaded triangles coloured by part (grey when not
/// selected, red where diagnosis found a defect, yellow for the chosen
/// record), and lines for the wireframe, defect edges and defect vertices.
fn geometry(
    session: &Session,
    wireframe: bool,
    defects: bool,
    record: Option<usize>,
) -> (Vec<Vertex>, Vec<Vertex>) {
    let scene = &session.scene;
    let model = &scene.model;
    let pos = |v: u32| {
        model
            .vertices
            .get(v as usize)
            .map_or([0.0; 3], |p| p.map(|c| c as f32))
    };
    let chosen: Option<&Record> = record.and_then(|i| scene.records.get(i));
    let mut defect_face = vec![false; model.triangles.len()];
    if defects {
        for &(f, _) in &scene.overlay.faces {
            if let Some(d) = defect_face.get_mut(f as usize) {
                *d = true;
            }
        }
    }
    let mut chosen_face = vec![false; model.triangles.len()];
    if let Some(r) = chosen {
        for &f in &r.location.faces {
            if let Some(c) = chosen_face.get_mut(f as usize) {
                *c = true;
            }
        }
    }
    let diag = {
        let [lo, hi] = scene.bounds;
        ((0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt() as f32).max(1e-6)
    };
    let mut tris = Vec::with_capacity(model.triangles.len() * 3);
    let mut lines = Vec::new();
    let line = |lines: &mut Vec<Vertex>, a: [f32; 3], b: [f32; 3], c: [f32; 4]| {
        for p in [a, b] {
            lines.push(Vertex {
                pos: p,
                normal: [0.0, 0.0, 1.0],
                color: c,
            });
        }
    };
    let edge_colour = [0.0, 0.0, 0.0, 0.35];
    for (t, tri) in model.triangles.iter().enumerate() {
        let Some(&part) = scene.triangle_part.get(t) else {
            continue;
        };
        if part == u32::MAX {
            continue;
        }
        let part = part as usize;
        let [a, b, c] = tri.vertices.map(pos);
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let normal = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let selected = session.selected.get(part).copied().unwrap_or(false);
        let colour = if chosen_face[t] {
            rgba(HIGHLIGHT, 1.0)
        } else if defect_face[t] {
            rgba(DEFECT, 1.0)
        } else if selected {
            rgba(part_colour(part), 1.0)
        } else {
            rgba(UNSELECTED, 1.0)
        };
        for p in [a, b, c] {
            tris.push(Vertex {
                pos: p,
                normal,
                color: colour,
            });
        }
        if wireframe {
            line(&mut lines, a, b, edge_colour);
            line(&mut lines, b, c, edge_colour);
            line(&mut lines, c, a, edge_colour);
        }
    }
    let marker = |lines: &mut Vec<Vertex>, p: [f32; 3], c: [f32; 4]| {
        let r = diag * 0.006;
        for k in 0..3 {
            let (mut lo, mut hi) = (p, p);
            lo[k] -= r;
            hi[k] += r;
            line(lines, lo, hi, c);
        }
    };
    if defects {
        let red = rgba(DEFECT, 1.0);
        for &([a, b], _) in &scene.overlay.edges {
            line(&mut lines, pos(a), pos(b), red);
        }
        for &(v, _) in &scene.overlay.vertices {
            marker(&mut lines, pos(v), red);
        }
    }
    if let Some(r) = chosen {
        let yellow = rgba(HIGHLIGHT, 1.0);
        for &v in &r.location.vertices {
            marker(&mut lines, pos(v), yellow);
        }
        if let Some(p) = r.location.point {
            marker(&mut lines, p.map(|c| c as f32), yellow);
        }
    }
    (tris, lines)
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.poll_loading();
        self.poll_run();
        if self.run.is_some() || self.loading.is_some() {
            ctx.request_repaint_after(Duration::from_millis(100));
        }
        // Files dropped on the window open like Open...
        let dropped: Option<PathBuf> = ctx.input(|i| {
            i.raw
                .dropped_files
                .first()
                .map(|f| f.path().to_path_buf())
                .filter(|p| !p.as_os_str().is_empty())
        });
        if let Some(path) = dropped {
            self.open(&path, &ctx);
        }

        let mut run_clicked = false;
        egui::Panel::top("toolbar")
            .frame(theme::bar_frame())
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    if ui.add(IconButton::new(Icon::Folder, "Open...")).clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("OBJ model or project", &["obj", "toml"])
                            .pick_file()
                        {
                            self.open(&path, &ctx);
                        }
                    }
                    ui.add_enabled_ui(self.session.is_some(), |ui| {
                        if ui
                            .add(IconButton::new(Icon::Save, "Save project..."))
                            .clicked()
                        {
                            if let Some(session) = &self.session {
                                if let Some(path) = rfd::FileDialog::new()
                                    .add_filter("Project", &["toml"])
                                    .set_file_name("project.toml")
                                    .save_file()
                                {
                                    self.status = match session.project().write(&path) {
                                        Ok(()) => format!("Saved {}", short_path(&path)),
                                        Err(e) => format!("Cannot save: {e}"),
                                    };
                                }
                            }
                        }
                        let can_run = self.run.is_none()
                            && self.session.as_ref().is_some_and(Session::any_selected);
                        if ui
                            .add_enabled(can_run, IconButton::primary(Icon::Play, "Run and export"))
                            .clicked()
                        {
                            run_clicked = true;
                        }
                    });
                    ui.separator();
                    if theme::pill(ui, &mut self.wireframe, "Wireframe").changed() {
                        self.dirty = true;
                    }
                    if theme::pill(ui, &mut self.defects, "Defects").changed() {
                        self.dirty = true;
                    }
                    theme::pill(ui, &mut self.clip, "Clip");
                    ui.add_enabled(
                        self.clip,
                        egui::Slider::new(&mut self.clip_at, 0.0..=1.0).show_value(false),
                    );
                });
            });
        egui::Panel::bottom("status")
            .frame(theme::bar_frame())
            .show(ui, |ui| {
                ui.label(theme::dim(self.status.as_str()));
            });
        egui::Panel::left("parts")
            .default_size(320.0)
            .frame(theme::side_frame())
            .show(ui, |ui| {
                self.parts_panel(ui);
            });
        egui::Panel::right("run")
            .default_size(400.0)
            .frame(theme::side_frame())
            .show(ui, |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("run panel")
                    .show(ui, |ui| self.run_panel(ui, &ctx));
            });
        egui::CentralPanel::default()
            .frame(egui::Frame::new())
            .show(ui, |ui| {
                self.viewport_panel(ui, frame);
            });
        if run_clicked {
            self.start_run(&ctx);
        }
    }
}

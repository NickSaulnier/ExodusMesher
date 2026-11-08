//! `mesher-gui`: pick parts of an OBJ file, run the worker, export. It holds
//! no meshing logic: models are read and diagnosed with the worker's own code
//! for display (on a background thread), and meshing runs in `mesher-worker`
//! through `mesher-client`, exactly as the CLI does. Every export saves a
//! project file that `mesher run project.toml` replays.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod render;

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

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("ExodusMesher")
            .with_inner_size([1280.0, 800.0]),
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

fn linear(c: [u8; 3], alpha: f32) -> [f32; 4] {
    let l = |v: u8| (f32::from(v) / 255.0).powf(2.2);
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
    run: Option<Run>,
    log: Vec<(Severity, String)>,
    status: String,
    worker: PathBuf,
    timeout_seconds: f64,
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>, start: Option<PathBuf>) -> App {
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
        self.status = format!("Reading {} ...", input.display());
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
                    scene.path.display(),
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
                self.status = format!("Reading {} failed", path.display());
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
            self.status = format!("Cannot create {}: {e}", session.output.display());
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
                        .map(|s| s.output.join(MESH_EXODUS_FILE).display().to_string())
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
        ui.heading("Parts");
        let Some(session) = &mut self.session else {
            ui.label("No model loaded.");
            return;
        };
        ui.horizontal(|ui| {
            if ui.button("All").clicked() {
                session.selected.iter_mut().for_each(|s| *s = true);
                self.dirty = true;
            }
            if ui.button("None").clicked() {
                session.selected.iter_mut().for_each(|s| *s = false);
                self.dirty = true;
            }
        });
        ui.separator();
        egui::ScrollArea::vertical().show(ui, |ui| {
            for i in 0..session.scene.parts.len() {
                let part = &session.scene.parts[i];
                let [r, g, b] = part_colour(i);
                let size = {
                    let [lo, hi] = part.bounds;
                    (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
                };
                let (tris, badges, name) = (part.triangles.len(), part.badges, part.name.clone());
                ui.horizontal(|ui| {
                    let mut on = session.selected[i];
                    if ui.checkbox(&mut on, "").changed() {
                        session.selected[i] = on;
                        self.dirty = true;
                    }
                    ui.label(RichText::new("■").color(Color32::from_rgb(r, g, b)));
                    ui.label(RichText::new(format!("{}. {name}", i + 1)).strong());
                });
                ui.horizontal(|ui| {
                    ui.add_space(24.0);
                    ui.label(format!("{tris} triangles, size {size:.3}"));
                    for (s, sev) in [Severity::Warning, Severity::Error, Severity::Fatal]
                        .iter()
                        .map(|&s| (badges[severity_slot(s)], s))
                    {
                        if s > 0 {
                            ui.label(
                                RichText::new(format!("{} {s}", short(sev)))
                                    .color(severity_colour(sev)),
                            );
                        }
                    }
                });
            }
        });
    }

    fn run_panel(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        ui.heading("Run");
        let running = self.run.is_some();
        let mut go = false;
        if let Some(session) = &mut self.session {
            let s = &mut session.settings;
            egui::Grid::new("settings").num_columns(2).show(ui, |ui| {
                ui.label("Element size").on_hover_text(
                    "Target element size, relative to each part's bounding-box diagonal",
                );
                ui.add(
                    egui::DragValue::new(&mut s.element_size)
                        .speed(0.001)
                        .range(1e-4..=1.0),
                );
                ui.end_row();
                ui.label("Max deviation").on_hover_text(
                    "Largest acceptable fallback deviation, relative to the diagonal",
                );
                ui.add(
                    egui::DragValue::new(&mut s.max_deviation)
                        .speed(0.001)
                        .range(1e-6..=1.0),
                );
                ui.end_row();
                ui.label("Strict");
                ui.checkbox(&mut s.strict, "conformal only, no fallback");
                ui.end_row();
                ui.label("No repair");
                ui.checkbox(&mut s.no_repair, "diagnose and report only");
                ui.end_row();
            });
            ui.horizontal(|ui| {
                ui.label("Output");
                if ui.button("Choose...").clicked() {
                    if let Some(dir) = rfd::FileDialog::new()
                        .set_directory(&session.output)
                        .pick_folder()
                    {
                        session.output = dir;
                    }
                }
            });
            ui.label(RichText::new(session.output.display().to_string()).small());
            ui.horizontal(|ui| {
                ui.add_enabled_ui(!running && session.any_selected(), |ui| {
                    if ui
                        .button("Run and export")
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
            if ui.button("Cancel").clicked() {
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
        ui.heading(format!("Records ({})", session.scene.records.len()));
        let mut clicked = None;
        egui::ScrollArea::vertical()
            .id_salt("records")
            .max_height(260.0)
            .show(ui, |ui| {
                for (i, r) in session.scene.records.iter().enumerate() {
                    let text = RichText::new(record_line(r))
                        .small()
                        .color(severity_colour(r.severity));
                    if ui.selectable_label(self.record == Some(i), text).clicked() {
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
        if let Some(r) = self.record.and_then(|i| session.scene.records.get(i)) {
            ui.separator();
            ui.label(
                RichText::new(r.code.as_str())
                    .strong()
                    .color(severity_colour(r.severity)),
            );
            ui.label(r.code.description());
            if let Some(part) = &r.part {
                ui.label(format!("Part: {part}"));
            }
            if let Some(hint) = &r.hint {
                ui.label(hint);
            }
            for (k, v) in r.measured.iter().chain(&r.threshold) {
                ui.label(RichText::new(format!("{k}: {v}")).small());
            }
            let lines = &r.location.obj_lines;
            if !lines.is_empty() {
                let shown: Vec<String> = lines.iter().take(20).map(u32::to_string).collect();
                ui.label(format!(
                    "OBJ lines: {}{}",
                    shown.join(", "),
                    if lines.len() > 20 { ", ..." } else { "" }
                ));
            }
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
                egui::FontId::proportional(16.0),
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
        let id = viewport.render(rs, size, &uniforms, [0.11, 0.12, 0.14, 1.0]);
        ui.painter().image(
            id,
            rect,
            egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
            Color32::WHITE,
        );
        ui.painter().text(
            rect.left_bottom() + egui::vec2(8.0, -8.0),
            egui::Align2::LEFT_BOTTOM,
            "drag: orbit   right-drag: pan   scroll: zoom   click: toggle part   double-click: fit",
            egui::FontId::proportional(12.0),
            Color32::from_gray(170),
        );
    }
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
    let edge_colour = [0.02, 0.02, 0.03, 0.55];
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
            linear(HIGHLIGHT, 1.0)
        } else if defect_face[t] {
            linear(DEFECT, 1.0)
        } else if selected {
            linear(part_colour(part), 1.0)
        } else {
            linear(UNSELECTED, 1.0)
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
        let red = linear(DEFECT, 1.0);
        for &([a, b], _) in &scene.overlay.edges {
            line(&mut lines, pos(a), pos(b), red);
        }
        for &(v, _) in &scene.overlay.vertices {
            marker(&mut lines, pos(v), red);
        }
    }
    if let Some(r) = chosen {
        let yellow = linear(HIGHLIGHT, 1.0);
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
        egui::Panel::top("toolbar").show(ui, |ui| {
            ui.horizontal(|ui| {
                if ui.button("Open...").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("OBJ model or project", &["obj", "toml"])
                        .pick_file()
                    {
                        self.open(&path, &ctx);
                    }
                }
                ui.add_enabled_ui(self.session.is_some(), |ui| {
                    if ui.button("Save project...").clicked() {
                        if let Some(session) = &self.session {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("Project", &["toml"])
                                .set_file_name("project.toml")
                                .save_file()
                            {
                                self.status = match session.project().write(&path) {
                                    Ok(()) => format!("Saved {}", path.display()),
                                    Err(e) => format!("Cannot save: {e}"),
                                };
                            }
                        }
                    }
                    let can_run = self.run.is_none()
                        && self.session.as_ref().is_some_and(Session::any_selected);
                    if ui
                        .add_enabled(can_run, egui::Button::new("Run and export"))
                        .clicked()
                    {
                        run_clicked = true;
                    }
                });
                ui.separator();
                if ui.checkbox(&mut self.wireframe, "Wireframe").changed() {
                    self.dirty = true;
                }
                if ui.checkbox(&mut self.defects, "Defects").changed() {
                    self.dirty = true;
                }
                ui.checkbox(&mut self.clip, "Clip");
                ui.add_enabled(
                    self.clip,
                    egui::Slider::new(&mut self.clip_at, 0.0..=1.0).show_value(false),
                );
            });
        });
        egui::Panel::bottom("status").show(ui, |ui| {
            ui.label(&self.status);
        });
        egui::Panel::left("parts")
            .default_size(280.0)
            .show(ui, |ui| {
                self.parts_panel(ui);
            });
        egui::Panel::right("run")
            .default_size(340.0)
            .show(ui, |ui| {
                self.run_panel(ui, &ctx);
            });
        egui::CentralPanel::default().show(ui, |ui| {
            self.viewport_panel(ui, frame);
        });
        if run_clicked {
            self.start_run(&ctx);
        }
    }
}

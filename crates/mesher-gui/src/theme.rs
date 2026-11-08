//! The window's look: larger text and controls, roomier spacing, rounded
//! corners and a dark slate palette with one accent colour.

use eframe::egui::{
    self, Color32, CornerRadius, FontFamily, FontId, Margin, Pos2, Rect, Response, Sense, Shape,
    Stroke, StrokeKind, TextStyle, Ui, Visuals, Widget, pos2, vec2,
};

/// The accent colour: primary buttons, selection, progress.
pub const ACCENT: Color32 = Color32::from_rgb(64, 132, 230);
const ACCENT_HOVER: Color32 = Color32::from_rgb(88, 152, 240);

const PANEL: Color32 = Color32::from_rgb(28, 31, 36);
const WINDOW: Color32 = Color32::from_rgb(34, 37, 43);
/// Background of inset areas (path fields, details).
pub const SUNKEN: Color32 = Color32::from_rgb(19, 21, 25);
/// Card fill for a selected row.
pub const CARD: Color32 = Color32::from_rgb(38, 42, 49);
/// Card fill under the pointer.
pub const CARD_HOVER: Color32 = Color32::from_rgb(47, 52, 61);
const FAINT: Color32 = Color32::from_rgb(36, 39, 45);
const BORDER: Color32 = Color32::from_rgb(50, 54, 62);
const BORDER_STRONG: Color32 = Color32::from_rgb(66, 71, 81);
const PILL_ON: Color32 = Color32::from_rgb(30, 54, 92);
const WIDGET: Color32 = Color32::from_rgb(45, 49, 57);
const WIDGET_HOVER: Color32 = Color32::from_rgb(58, 63, 73);
const WIDGET_ACTIVE: Color32 = Color32::from_rgb(68, 74, 86);
const TEXT: Color32 = Color32::from_rgb(222, 226, 232);
pub const TEXT_DIM: Color32 = Color32::from_rgb(160, 166, 176);

const RADIUS: u8 = 6;

/// Install the theme on a context.
pub fn apply(ctx: &egui::Context) {
    ctx.set_theme(egui::Theme::Dark);
    ctx.style_mut_of(egui::Theme::Dark, |style| {
        style.text_styles = [
            (
                TextStyle::Small,
                FontId::new(12.5, FontFamily::Proportional),
            ),
            (TextStyle::Body, FontId::new(15.5, FontFamily::Proportional)),
            (
                TextStyle::Button,
                FontId::new(15.5, FontFamily::Proportional),
            ),
            (
                TextStyle::Heading,
                FontId::new(21.0, FontFamily::Proportional),
            ),
            (
                TextStyle::Monospace,
                FontId::new(14.0, FontFamily::Monospace),
            ),
        ]
        .into();

        let s = &mut style.spacing;
        s.item_spacing = vec2(10.0, 8.0);
        s.button_padding = vec2(14.0, 7.0);
        s.interact_size = vec2(44.0, 32.0);
        s.window_margin = Margin::same(14);
        s.menu_margin = Margin::same(8);
        s.indent = 22.0;
        s.icon_width = 20.0;
        s.icon_width_inner = 11.0;
        s.icon_spacing = 8.0;
        s.slider_width = 160.0;
        s.combo_width = 160.0;
        s.text_edit_width = 220.0;
        s.scroll.bar_width = 10.0;

        style.visuals = visuals();
    });
}

fn visuals() -> Visuals {
    let mut v = Visuals::dark();
    v.panel_fill = PANEL;
    v.window_fill = WINDOW;
    v.extreme_bg_color = SUNKEN;
    v.faint_bg_color = FAINT;
    v.code_bg_color = SUNKEN;
    v.window_stroke = Stroke::new(1.0, BORDER);
    v.window_corner_radius = CornerRadius::same(10);
    v.menu_corner_radius = CornerRadius::same(8);
    v.hyperlink_color = ACCENT_HOVER;
    v.selection.bg_fill = Color32::from_rgb(44, 92, 165);
    v.selection.stroke = Stroke::new(1.0, Color32::WHITE);

    let w = &mut v.widgets;
    w.noninteractive.bg_fill = PANEL;
    w.noninteractive.weak_bg_fill = PANEL;
    w.noninteractive.bg_stroke = Stroke::new(1.0, BORDER);
    w.noninteractive.fg_stroke = Stroke::new(1.0, TEXT);
    w.noninteractive.corner_radius = CornerRadius::same(RADIUS);

    for (state, fill, text) in [
        (&mut w.inactive, WIDGET, TEXT_DIM),
        (&mut w.hovered, WIDGET_HOVER, TEXT),
        (&mut w.active, WIDGET_ACTIVE, Color32::WHITE),
        (&mut w.open, WIDGET_HOVER, TEXT),
    ] {
        state.bg_fill = fill;
        state.weak_bg_fill = fill;
        state.bg_stroke = Stroke::NONE;
        state.fg_stroke = Stroke::new(1.5, text);
        state.corner_radius = CornerRadius::same(RADIUS);
        state.expansion = 0.0;
    }
    w.inactive.fg_stroke = Stroke::new(1.0, TEXT);
    w.hovered.bg_stroke = Stroke::new(1.0, Color32::from_rgb(80, 86, 98));
    w.active.bg_stroke = Stroke::new(1.0, ACCENT);
    v
}

/// Frame for the top toolbar and the status bar.
pub fn bar_frame() -> egui::Frame {
    egui::Frame::new()
        .inner_margin(Margin::symmetric(14, 10))
        .fill(PANEL)
}

/// Frame for the side panels.
pub fn side_frame() -> egui::Frame {
    egui::Frame::new()
        .inner_margin(Margin::same(16))
        .fill(PANEL)
}

/// Dimmed secondary text.
#[must_use]
pub fn dim(text: impl Into<String>) -> egui::RichText {
    egui::RichText::new(text).color(TEXT_DIM)
}

/// Small vector icons, drawn in a 16 × 16 box so they need no icon font.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Icon {
    Folder,
    Save,
    Play,
    Stop,
    Check,
    Cross,
}

impl Icon {
    /// Paint the icon into `rect` (square) in `color`.
    pub fn paint(self, painter: &egui::Painter, rect: Rect, color: Color32) {
        let k = rect.width() / 16.0;
        let p = |x: f32, y: f32| -> Pos2 { rect.min + vec2(x * k, y * k) };
        let stroke = Stroke::new(1.6 * k.max(1.0), color);
        match self {
            Icon::Folder => {
                let outline = vec![
                    p(1.5, 3.5),
                    p(6.0, 3.5),
                    p(7.5, 5.0),
                    p(14.5, 5.0),
                    p(14.5, 13.0),
                    p(1.5, 13.0),
                ];
                painter.add(Shape::closed_line(outline, stroke));
                painter.line_segment([p(1.5, 7.0), p(14.5, 7.0)], stroke);
            }
            Icon::Save => {
                painter.rect_stroke(
                    Rect::from_min_max(p(2.0, 2.0), p(14.0, 14.0)),
                    2.0 * k,
                    stroke,
                    StrokeKind::Middle,
                );
                painter.rect_stroke(
                    Rect::from_min_max(p(5.0, 2.0), p(11.0, 6.0)),
                    0.0,
                    stroke,
                    StrokeKind::Middle,
                );
                painter.rect_filled(Rect::from_min_max(p(5.0, 9.5), p(11.0, 14.0)), 1.0, color);
            }
            Icon::Play => {
                painter.add(Shape::convex_polygon(
                    vec![p(4.0, 2.5), p(13.5, 8.0), p(4.0, 13.5)],
                    color,
                    Stroke::NONE,
                ));
            }
            Icon::Stop => {
                painter.rect_filled(
                    Rect::from_min_max(p(3.5, 3.5), p(12.5, 12.5)),
                    2.0 * k,
                    color,
                );
            }
            Icon::Check => {
                painter.add(Shape::line(
                    vec![p(3.0, 8.5), p(6.5, 12.0), p(13.0, 4.5)],
                    stroke,
                ));
            }
            Icon::Cross => {
                painter.line_segment([p(4.0, 4.0), p(12.0, 12.0)], stroke);
                painter.line_segment([p(12.0, 4.0), p(4.0, 12.0)], stroke);
            }
        }
    }
}

/// A button with an icon before its text; `primary` fills it with the
/// accent colour.
pub struct IconButton<'a> {
    icon: Icon,
    text: &'a str,
    primary: bool,
}

impl<'a> IconButton<'a> {
    pub fn new(icon: Icon, text: &'a str) -> Self {
        IconButton {
            icon,
            text,
            primary: false,
        }
    }

    pub fn primary(icon: Icon, text: &'a str) -> Self {
        IconButton {
            icon,
            text,
            primary: true,
        }
    }
}

impl Widget for IconButton<'_> {
    fn ui(self, ui: &mut Ui) -> Response {
        const ICON: f32 = 16.0;
        let font = TextStyle::Button.resolve(ui.style());
        let galley = ui
            .painter()
            .layout_no_wrap(self.text.to_owned(), font, Color32::PLACEHOLDER);
        let pad = ui.spacing().button_padding;
        let gap = if self.text.is_empty() { 0.0 } else { 8.0 };
        let size = vec2(
            pad.x * 2.0 + ICON + gap + galley.size().x,
            (galley.size().y + pad.y * 2.0).max(ui.spacing().interact_size.y),
        );
        let (rect, response) = ui.allocate_exact_size(size, Sense::click());
        if ui.is_rect_visible(rect) {
            let visuals = ui.style().interact(&response);
            let (fill, fg, stroke) = if self.primary {
                let fill = if response.hovered() || response.is_pointer_button_down_on() {
                    ACCENT_HOVER
                } else {
                    ACCENT
                };
                (fill, Color32::WHITE, Stroke::NONE)
            } else {
                (
                    visuals.weak_bg_fill,
                    visuals.text_color(),
                    visuals.bg_stroke,
                )
            };
            let painter = ui.painter();
            painter.rect(
                rect,
                visuals.corner_radius,
                fill,
                stroke,
                StrokeKind::Inside,
            );
            let icon = Rect::from_min_size(
                pos2(rect.left() + pad.x, rect.center().y - ICON * 0.5),
                vec2(ICON, ICON),
            );
            self.icon.paint(painter, icon, fg);
            let text = pos2(icon.right() + gap, rect.center().y - galley.size().y * 0.5);
            painter.galley(text, galley, fg);
        }
        response
    }
}

/// A rounded on/off chip with a status dot, for view options.
pub fn pill(ui: &mut Ui, on: &mut bool, text: &str) -> Response {
    const DOT: f32 = 8.0;
    let font = TextStyle::Button.resolve(ui.style());
    let galley = ui
        .painter()
        .layout_no_wrap(text.to_owned(), font, Color32::PLACEHOLDER);
    let size = vec2(
        14.0 * 2.0 + DOT + 8.0 + galley.size().x,
        ui.spacing().interact_size.y,
    );
    let (rect, mut response) = ui.allocate_exact_size(size, Sense::click());
    if response.clicked() {
        *on = !*on;
        response.mark_changed();
    }
    if ui.is_rect_visible(rect) {
        let (fill, stroke, fg) = if *on {
            (PILL_ON, Stroke::new(1.0, ACCENT), Color32::WHITE)
        } else if response.hovered() {
            (WIDGET_HOVER, Stroke::new(1.0, BORDER_STRONG), TEXT)
        } else {
            (
                Color32::TRANSPARENT,
                Stroke::new(1.0, BORDER_STRONG),
                TEXT_DIM,
            )
        };
        let painter = ui.painter();
        let round = CornerRadius::same((rect.height() * 0.5) as u8);
        painter.rect(rect, round, fill, stroke, StrokeKind::Inside);
        let dot = pos2(rect.left() + 14.0 + DOT * 0.5, rect.center().y);
        if *on {
            painter.circle_filled(dot, DOT * 0.5, ACCENT_HOVER);
        } else {
            painter.circle_stroke(dot, DOT * 0.5 - 0.5, Stroke::new(1.2, TEXT_DIM));
        }
        let text = pos2(
            rect.left() + 14.0 + DOT + 8.0,
            rect.center().y - galley.size().y * 0.5,
        );
        painter.galley(text, galley, fg);
    }
    response
}

/// A sliding on/off switch.
pub fn switch(ui: &mut Ui, on: &mut bool) -> Response {
    let size = vec2(38.0, 22.0);
    let (rect, mut response) = ui.allocate_exact_size(size, Sense::click());
    if response.clicked() {
        *on = !*on;
        response.mark_changed();
    }
    if ui.is_rect_visible(rect) {
        let t = ui.ctx().animate_bool_responsive(response.id, *on);
        let fill = if *on {
            ACCENT
        } else if response.hovered() {
            WIDGET_ACTIVE
        } else {
            WIDGET_HOVER
        };
        let painter = ui.painter();
        painter.rect_filled(rect, CornerRadius::same(11), fill);
        let r = rect.height() * 0.5 - 3.0;
        let x = egui::lerp((rect.left() + 3.0 + r)..=(rect.right() - 3.0 - r), t);
        painter.circle_filled(pos2(x, rect.center().y), r, Color32::WHITE);
    }
    response
}

/// A small tinted label, such as a severity tag. With `active` false it is
/// drawn dimmed; it senses clicks so it can act as a filter.
pub fn chip(ui: &mut Ui, text: &str, color: Color32, active: bool) -> Response {
    let font = TextStyle::Small.resolve(ui.style());
    let galley = ui
        .painter()
        .layout_no_wrap(text.to_owned(), font, Color32::PLACEHOLDER);
    let size = vec2((galley.size().x + 14.0).max(48.0), galley.size().y + 6.0);
    let (rect, response) = ui.allocate_exact_size(size, Sense::click());
    if ui.is_rect_visible(rect) {
        let (fill, fg) = if active {
            (color.gamma_multiply(0.22), color)
        } else {
            (Color32::TRANSPARENT, TEXT_DIM)
        };
        let painter = ui.painter();
        painter.rect(
            rect,
            CornerRadius::same(5),
            fill,
            Stroke::new(
                1.0,
                if active {
                    color.gamma_multiply(0.5)
                } else {
                    BORDER_STRONG
                },
            ),
            StrokeKind::Inside,
        );
        let text = rect.center() - galley.size() * 0.5;
        painter.galley(text, galley, fg);
    }
    response
}

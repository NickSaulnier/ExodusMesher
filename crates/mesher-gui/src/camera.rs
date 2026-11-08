//! An orbit camera: a target point, a distance, and yaw and pitch angles.
//! Matrices are column-major (as WGSL expects) with a [0, 1] depth range.

pub type Vec3 = [f32; 3];
pub type Mat4 = [[f32; 4]; 4];

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: Vec3, s: f32) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: Vec3, b: Vec3) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn normalize(a: Vec3) -> Vec3 {
    let n = dot(a, a).sqrt();
    if n > 0.0 {
        scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 1.0]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    pub target: Vec3,
    pub distance: f32,
    /// Rotation about +Z, radians.
    pub yaw: f32,
    /// Elevation above the XY plane, radians (kept inside ±89°).
    pub pitch: f32,
    /// Vertical field of view, radians.
    pub fov_y: f32,
}

const PITCH_LIMIT: f32 = 1.55;

impl Camera {
    /// Frame a bounding box `[min, max]`.
    #[must_use]
    pub fn fit(bounds: [[f64; 3]; 2]) -> Camera {
        let [lo, hi] = bounds.map(|p| p.map(|c| c as f32));
        let target = scale(add(lo, hi), 0.5);
        let radius = (dot(sub(hi, lo), sub(hi, lo)).sqrt() * 0.5).max(1e-6);
        let fov_y = 45f32.to_radians();
        Camera {
            target,
            distance: radius / (fov_y * 0.5).sin() * 1.1,
            yaw: -0.6,
            pitch: 0.45,
            fov_y,
        }
    }

    /// Unit vector from the target towards the eye.
    #[must_use]
    pub fn back(&self) -> Vec3 {
        let (cp, sp) = (self.pitch.cos(), self.pitch.sin());
        [cp * self.yaw.cos(), cp * self.yaw.sin(), sp]
    }

    #[must_use]
    pub fn eye(&self) -> Vec3 {
        add(self.target, scale(self.back(), self.distance))
    }

    /// Camera right and up vectors.
    fn basis(&self) -> (Vec3, Vec3, Vec3) {
        let f = scale(self.back(), -1.0);
        let right = normalize(cross(f, [0.0, 0.0, 1.0]));
        let up = cross(right, f);
        (f, right, up)
    }

    /// Rotate by a mouse drag of `dx`, `dy` pixels.
    pub fn orbit(&mut self, dx: f32, dy: f32) {
        self.yaw -= dx * 0.01;
        self.pitch = (self.pitch + dy * 0.01).clamp(-PITCH_LIMIT, PITCH_LIMIT);
    }

    /// Slide the target by a drag of `dx`, `dy` pixels in a viewport
    /// `height` pixels tall, so the point under the cursor follows it.
    pub fn pan(&mut self, dx: f32, dy: f32, height: f32) {
        let (_, right, up) = self.basis();
        let per_pixel = 2.0 * self.distance * (self.fov_y * 0.5).tan() / height.max(1.0);
        self.target = add(
            self.target,
            add(scale(right, -dx * per_pixel), scale(up, dy * per_pixel)),
        );
    }

    /// Move towards the target by `factor` (< 1 zooms in).
    pub fn zoom(&mut self, factor: f32) {
        self.distance = (self.distance * factor).max(1e-6);
    }

    /// Near and far planes around the target.
    fn depth_range(&self) -> (f32, f32) {
        (self.distance * 0.01, self.distance * 100.0)
    }

    /// Projection × view, column-major, depth in [0, 1].
    #[must_use]
    pub fn view_proj(&self, aspect: f32) -> Mat4 {
        let (f, right, up) = self.basis();
        let eye = self.eye();
        // View: rows are right, up, -forward.
        let view = [
            [right[0], up[0], -f[0], 0.0],
            [right[1], up[1], -f[1], 0.0],
            [right[2], up[2], -f[2], 0.0],
            [-dot(right, eye), -dot(up, eye), dot(f, eye), 1.0],
        ];
        let (near, far) = self.depth_range();
        let t = 1.0 / (self.fov_y * 0.5).tan();
        let proj = [
            [t / aspect.max(1e-6), 0.0, 0.0, 0.0],
            [0.0, t, 0.0, 0.0],
            [0.0, 0.0, far / (near - far), -1.0],
            [0.0, 0.0, near * far / (near - far), 0.0],
        ];
        mul(proj, view)
    }

    /// The ray through a point of the viewport, in normalized device
    /// coordinates (x right, y up, both in [-1, 1]): origin and direction.
    #[must_use]
    pub fn ray(&self, x: f32, y: f32, aspect: f32) -> ([f64; 3], [f64; 3]) {
        let (f, right, up) = self.basis();
        let t = (self.fov_y * 0.5).tan();
        let dir = normalize(add(f, add(scale(right, x * t * aspect), scale(up, y * t))));
        (self.eye().map(f64::from), dir.map(f64::from))
    }
}

/// `a × b` for column-major matrices.
fn mul(a: Mat4, b: Mat4) -> Mat4 {
    let mut out = [[0.0; 4]; 4];
    for (c, col) in out.iter_mut().enumerate() {
        for (r, v) in col.iter_mut().enumerate() {
            *v = (0..4).map(|k| a[k][r] * b[c][k]).sum();
        }
    }
    out
}

/// Apply a column-major matrix to a point; returns clip coordinates.
#[must_use]
pub fn transform(m: Mat4, p: Vec3) -> [f32; 4] {
    let v = [p[0], p[1], p[2], 1.0];
    let mut out = [0.0; 4];
    for (r, o) in out.iter_mut().enumerate() {
        *o = (0..4).map(|c| m[c][r] * v[c]).sum();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn target_projects_to_the_centre_and_inside_the_depth_range() {
        let cam = Camera::fit([[0.0, 0.0, 0.0], [2.0, 1.0, 1.0]]);
        let c = transform(cam.view_proj(1.5), cam.target);
        assert!(c[3] > 0.0);
        let (x, y, z) = (c[0] / c[3], c[1] / c[3], c[2] / c[3]);
        assert!(x.abs() < 1e-5 && y.abs() < 1e-5, "{x} {y}");
        assert!((0.0..1.0).contains(&z), "{z}");
    }

    #[test]
    fn rays_pass_through_what_projects_there() {
        let mut cam = Camera::fit([[0.0; 3], [1.0; 3]]);
        cam.orbit(40.0, -25.0);
        cam.pan(10.0, 5.0, 600.0);
        let aspect = 1.3;
        let p = [0.8f32, 0.1, 0.3];
        let c = transform(cam.view_proj(aspect), p);
        let (x, y) = (c[0] / c[3], c[1] / c[3]);
        let (o, d) = cam.ray(x, y, aspect);
        // Distance from p to the ray is tiny.
        let v = [
            f64::from(p[0]) - o[0],
            f64::from(p[1]) - o[1],
            f64::from(p[2]) - o[2],
        ];
        let t = v[0] * d[0] + v[1] * d[1] + v[2] * d[2];
        let off = [v[0] - t * d[0], v[1] - t * d[1], v[2] - t * d[2]];
        let dist = (off[0] * off[0] + off[1] * off[1] + off[2] * off[2]).sqrt();
        assert!(t > 0.0 && dist < 1e-4, "{dist}");
    }

    #[test]
    fn pitch_is_clamped_and_zoom_stays_positive() {
        let mut cam = Camera::fit([[0.0; 3], [1.0; 3]]);
        cam.orbit(0.0, 1e6);
        assert!(cam.pitch <= PITCH_LIMIT);
        cam.zoom(0.0);
        assert!(cam.distance > 0.0);
    }
}

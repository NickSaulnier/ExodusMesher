//! Clean shapes and one-defect injectors for the synthetic suite.
//!
//! Shapes are triangle meshes with shared vertices, outward-oriented and free
//! of defects. A [`Scene`] writes one or more of them as OBJ text, optionally
//! under `o` names; triangle `t` of the scene is the `t`-th `f` line, so test
//! expectations can name faces by index.

use std::collections::BTreeMap;

pub type Point = [f64; 3];

#[derive(Debug, Clone, PartialEq)]
pub struct Solid {
    pub points: Vec<Point>,
    pub tris: Vec<[u32; 3]>,
}

impl Solid {
    /// Axis-aligned box `[min, max]`, 12 outward triangles.
    #[must_use]
    pub fn cuboid(min: Point, max: Point) -> Solid {
        let p = |x: usize, y: usize, z: usize| {
            [
                if x == 0 { min[0] } else { max[0] },
                if y == 0 { min[1] } else { max[1] },
                if z == 0 { min[2] } else { max[2] },
            ]
        };
        let points = vec![
            p(0, 0, 0),
            p(1, 0, 0),
            p(1, 1, 0),
            p(0, 1, 0),
            p(0, 0, 1),
            p(1, 0, 1),
            p(1, 1, 1),
            p(0, 1, 1),
        ];
        let quads = [
            [0, 3, 2, 1], // z = min
            [4, 5, 6, 7], // z = max
            [0, 1, 5, 4], // y = min
            [2, 3, 7, 6], // y = max
            [0, 4, 7, 3], // x = min
            [1, 2, 6, 5], // x = max
        ];
        let tris = quads
            .iter()
            .flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]])
            .collect();
        Solid { points, tris }
    }

    #[must_use]
    pub fn unit_cube() -> Solid {
        Solid::cuboid([0.0; 3], [1.0; 3])
    }

    /// The prism of a simple, counter-clockwise polygon between `z0` and
    /// `z1` (outward). Its caps are fanned from vertex `fan`, which the
    /// caller picks so that every fan triangle lies inside the polygon.
    #[must_use]
    pub fn prism(polygon: &[[f64; 2]], z0: f64, z1: f64, fan: usize) -> Solid {
        let n = polygon.len() as u32;
        let points = [z0, z1]
            .iter()
            .flat_map(|&z| polygon.iter().map(move |&[x, y]| [x, y, z]))
            .collect();
        let at = |k: u32| (fan as u32 + k) % n;
        let mut tris = Vec::new();
        for k in 1..n.saturating_sub(1) {
            let (a, b, c) = (at(0), at(k), at(k + 1));
            tris.push([a, c, b]);
            tris.push([n + a, n + b, n + c]);
        }
        for i in 0..n {
            let j = (i + 1) % n;
            tris.push([i, j, n + j]);
            tris.push([i, n + j, n + i]);
        }
        Solid { points, tris }
    }

    /// Surface of `[0, n]^3` from unit quads, sharing vertices along edges.
    #[must_use]
    pub fn grid_cube(n: i64) -> Solid {
        let mut ids: BTreeMap<[i64; 3], u32> = BTreeMap::new();
        let mut points = Vec::new();
        let mut tris = Vec::new();
        let (x, y, z) = ([1, 0, 0], [0, 1, 0], [0, 0, 1]);
        let faces = [
            ([0, 0, 0], y, x),
            ([0, 0, n], x, y),
            ([0, 0, 0], x, z),
            ([0, n, 0], z, x),
            ([0, 0, 0], z, y),
            ([n, 0, 0], y, z),
        ];
        for (o, u, v) in faces {
            let mut id = |i: i64, j: i64| {
                let p = [0, 1, 2].map(|k| o[k] + i * u[k] + j * v[k]);
                *ids.entry(p).or_insert_with(|| {
                    points.push(p.map(|c| c as f64));
                    (points.len() - 1) as u32
                })
            };
            for i in 0..n {
                for j in 0..n {
                    let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
                    tris.push([a, b, c]);
                    tris.push([a, c, d]);
                }
            }
        }
        Solid { points, tris }
    }

    /// Unit icosphere.
    #[must_use]
    pub fn icosphere(subdivisions: u32) -> Solid {
        let t = (1.0 + 5f64.sqrt()) / 2.0;
        let mut points: Vec<Point> = vec![
            [-1.0, t, 0.0],
            [1.0, t, 0.0],
            [-1.0, -t, 0.0],
            [1.0, -t, 0.0],
            [0.0, -1.0, t],
            [0.0, 1.0, t],
            [0.0, -1.0, -t],
            [0.0, 1.0, -t],
            [t, 0.0, -1.0],
            [t, 0.0, 1.0],
            [-t, 0.0, -1.0],
            [-t, 0.0, 1.0],
        ];
        let mut tris: Vec<[u32; 3]> = vec![
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];
        for _ in 0..subdivisions {
            let mut mid: BTreeMap<(u32, u32), u32> = BTreeMap::new();
            let mut next = Vec::with_capacity(tris.len() * 4);
            for [a, b, c] in tris {
                let mut m = |i: u32, j: u32| {
                    *mid.entry((i.min(j), i.max(j))).or_insert_with(|| {
                        let (p, q) = (points[i as usize], points[j as usize]);
                        points.push([0, 1, 2].map(|k| (p[k] + q[k]) / 2.0));
                        (points.len() - 1) as u32
                    })
                };
                let (ab, bc, ca) = (m(a, b), m(b, c), m(c, a));
                next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
            }
            tris = next;
        }
        for p in &mut points {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            *p = p.map(|c| c / r);
        }
        Solid { points, tris }
    }

    /// Torus with major radius `big`, minor radius `small`. The angles use
    /// the pure-Rust `libm`, not the platform's maths library, whose last
    /// bits differ between Windows and Linux: the golden reports must see the
    /// same input everywhere.
    #[must_use]
    pub fn torus(big: f64, small: f64, nu: u32, nv: u32) -> Solid {
        let mut points = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let u = f64::from(i) / f64::from(nu) * std::f64::consts::TAU;
                let v = f64::from(j) / f64::from(nv) * std::f64::consts::TAU;
                let a = big + small * libm::cos(v);
                points.push([a * libm::cos(u), a * libm::sin(u), small * libm::sin(v)]);
            }
        }
        let id = |i: u32, j: u32| (i % nu) * nv + (j % nv);
        let mut tris = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            }
        }
        Solid { points, tris }
    }

    /// Reverse every triangle.
    #[must_use]
    pub fn flipped(mut self) -> Solid {
        for t in &mut self.tris {
            t.swap(1, 2);
        }
        self
    }

    /// Join two solids into one (vertex lists concatenated, then exactly
    /// coincident points shared so the join itself adds no duplicate vertices).
    #[must_use]
    pub fn union(mut self, other: &Solid) -> Solid {
        let offset = self.points.len() as u32;
        self.points.extend(&other.points);
        self.tris
            .extend(other.tris.iter().map(|t| t.map(|v| v + offset)));
        self.share_coincident()
    }

    /// Replace exactly coincident points by the first of them.
    #[must_use]
    pub fn share_coincident(self) -> Solid {
        let mut first: BTreeMap<[u64; 3], u32> = BTreeMap::new();
        let mut remap = Vec::with_capacity(self.points.len());
        let mut points = Vec::new();
        for p in &self.points {
            let key = p.map(|c| (c + 0.0).to_bits());
            let id = *first.entry(key).or_insert_with(|| {
                points.push(*p);
                (points.len() - 1) as u32
            });
            remap.push(id);
        }
        let tris = self
            .tris
            .iter()
            .map(|t| t.map(|v| remap[v as usize]))
            .collect();
        Solid { points, tris }
    }

    /// Index of the point at exactly `p`.
    #[must_use]
    pub fn find(&self, p: Point) -> Option<u32> {
        self.points.iter().position(|&q| q == p).map(|i| i as u32)
    }

    /// Bounding-box diagonal.
    #[must_use]
    pub fn diagonal(&self) -> f64 {
        let lo = [0, 1, 2].map(|k| {
            self.points
                .iter()
                .map(|p| p[k])
                .fold(f64::INFINITY, f64::min)
        });
        let hi = [0, 1, 2].map(|k| {
            self.points
                .iter()
                .map(|p| p[k])
                .fold(f64::NEG_INFINITY, f64::max)
        });
        (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
    }
}

/// Solids written into one OBJ file, each optionally under an `o` name.
#[derive(Debug, Clone, Default)]
pub struct Scene {
    pub solids: Vec<(Option<String>, Solid)>,
}

impl Scene {
    #[must_use]
    pub fn single(solid: Solid) -> Scene {
        Scene {
            solids: vec![(None, solid)],
        }
    }

    #[must_use]
    pub fn with(mut self, name: Option<&str>, solid: Solid) -> Scene {
        self.solids.push((name.map(str::to_owned), solid));
        self
    }

    /// Apply `p -> scale * p + offset` to every point.
    #[must_use]
    pub fn transformed(&self, scale: f64, offset: f64) -> Scene {
        Scene {
            solids: self
                .solids
                .iter()
                .map(|(n, s)| {
                    let mut s = s.clone();
                    for p in &mut s.points {
                        *p = p.map(|c| c * scale + offset);
                    }
                    (n.clone(), s)
                })
                .collect(),
        }
    }

    /// OBJ text: all vertices of a solid, then its faces, solid by solid.
    #[must_use]
    pub fn to_obj(&self) -> String {
        let mut out = String::from("# synthetic\n");
        let mut base = 1u32;
        for (name, solid) in &self.solids {
            if let Some(name) = name {
                out.push_str(&format!("o {name}\n"));
            }
            for p in &solid.points {
                out.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
            }
            for t in &solid.tris {
                out.push_str(&format!(
                    "f {} {} {}\n",
                    t[0] + base,
                    t[1] + base,
                    t[2] + base
                ));
            }
            base += solid.points.len() as u32;
        }
        out
    }

    /// Scene-wide index of solid `s`'s first vertex and first triangle.
    #[must_use]
    pub fn offsets(&self, s: usize) -> (u32, u32) {
        self.solids[..s].iter().fold((0, 0), |(v, t), (_, solid)| {
            (v + solid.points.len() as u32, t + solid.tris.len() as u32)
        })
    }
}

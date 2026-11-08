//! Helpers for oracle comparisons. TetGen itself is only a dev-dependency.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use mesher_geom::Point3;

/// A triangle surface with shared vertices.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Surface {
    pub points: Vec<Point3>,
    pub triangles: Vec<[u32; 3]>,
}

/// Read an STL file (binary or ASCII) and weld vertices whose coordinates are
/// bitwise equal (STL repeats every vertex per triangle). No other repair:
/// near-coincident vertices stay distinct, so the mesher's diagnosis sees the
/// geometry as published.
///
/// # Errors
/// A description of the first problem when the bytes are not a readable STL.
pub fn read_stl(bytes: &[u8]) -> Result<Surface, String> {
    let corners = if is_binary_stl(bytes) {
        binary_corners(bytes)?
    } else {
        ascii_corners(bytes)?
    };
    let mut index: BTreeMap<[u64; 3], u32> = BTreeMap::new();
    let mut surface = Surface::default();
    for tri in corners.chunks_exact(3) {
        let mut t = [0u32; 3];
        for (k, p) in tri.iter().enumerate() {
            // +0.0 and -0.0 are the same point.
            let key = p.map(|c| if c == 0.0 { 0 } else { c.to_bits() });
            t[k] = *index.entry(key).or_insert_with(|| {
                surface.points.push(*p);
                (surface.points.len() - 1) as u32
            });
        }
        surface.triangles.push(t);
    }
    Ok(surface)
}

/// Binary STL: an 80-byte header, a triangle count, then 50 bytes per
/// triangle. Some binary files start with "solid" too, so the size decides.
fn is_binary_stl(bytes: &[u8]) -> bool {
    match bytes.get(80..84) {
        Some(n) => {
            let n = u32::from_le_bytes([n[0], n[1], n[2], n[3]]) as usize;
            n.checked_mul(50).and_then(|b| b.checked_add(84)) == Some(bytes.len())
        }
        None => false,
    }
}

fn binary_corners(bytes: &[u8]) -> Result<Vec<Point3>, String> {
    let mut out = Vec::new();
    for record in bytes.get(84..).unwrap_or(&[]).chunks_exact(50) {
        // 12 bytes of normal, then three vertices of three f32 each.
        for v in 0..3 {
            let mut p = [0.0; 3];
            for (k, c) in p.iter_mut().enumerate() {
                let at = 12 + 12 * v + 4 * k;
                let b = record
                    .get(at..at + 4)
                    .ok_or_else(|| "truncated binary STL record".to_owned())?;
                *c = f64::from(f32::from_le_bytes([b[0], b[1], b[2], b[3]]));
            }
            if p.iter().any(|c| !c.is_finite()) {
                return Err("non-finite coordinate".to_owned());
            }
            out.push(p);
        }
    }
    Ok(out)
}

fn ascii_corners(bytes: &[u8]) -> Result<Vec<Point3>, String> {
    let text = std::str::from_utf8(bytes).map_err(|e| format!("ASCII STL is not UTF-8: {e}"))?;
    let mut out = Vec::new();
    let mut words = text.split_ascii_whitespace();
    while let Some(w) = words.next() {
        if !w.eq_ignore_ascii_case("vertex") {
            continue;
        }
        let mut p = [0.0; 3];
        for c in &mut p {
            *c = words
                .next()
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|c| c.is_finite())
                .ok_or_else(|| "bad vertex coordinate".to_owned())?;
        }
        out.push(p);
    }
    if out.len() % 3 != 0 {
        return Err(format!(
            "{} vertices is not a whole number of facets",
            out.len()
        ));
    }
    Ok(out)
}

/// The surface as OBJ text, with a single object named `name`. Coordinates
/// are written with Rust's shortest round-trip formatting, so they read back
/// exactly.
#[must_use]
pub fn to_obj(surface: &Surface, name: &str) -> String {
    let mut s = String::with_capacity(32 * (surface.points.len() + surface.triangles.len()));
    let _ = writeln!(s, "o {name}");
    for p in &surface.points {
        let _ = writeln!(s, "v {:?} {:?} {:?}", p[0], p[1], p[2]);
    }
    for t in &surface.triangles {
        let _ = writeln!(s, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1);
    }
    s
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    const TET: &str = "solid t
facet normal 0 0 -1
 outer loop
  vertex 0 0 0
  vertex 0 1 0
  vertex 1 0 0
 endloop
endfacet
facet normal 0 -1 0
 outer loop
  vertex 0 0 0
  vertex 1 0 0
  vertex 0 0 1
 endloop
endfacet
endsolid t
";

    #[test]
    fn ascii_stl_is_welded() {
        let s = read_stl(TET.as_bytes()).unwrap();
        assert_eq!(s.points.len(), 4);
        assert_eq!(s.triangles, vec![[0, 1, 2], [0, 2, 3]]);
    }

    #[test]
    fn binary_stl_matches_ascii() {
        let ascii = read_stl(TET.as_bytes()).unwrap();
        let mut b = b"solid but actually binary".to_vec();
        b.resize(80, 0);
        b.extend(2u32.to_le_bytes());
        for t in &ascii.triangles {
            b.extend([0u8; 12]);
            for &v in t {
                for c in ascii.points[v as usize] {
                    b.extend((c as f32).to_le_bytes());
                }
            }
            b.extend([0u8; 2]);
        }
        assert_eq!(read_stl(&b).unwrap(), ascii);
        let obj = to_obj(&ascii, "t");
        assert!(obj.contains("v 1.0 0.0 0.0") && obj.contains("f 1 3 4"));
    }
}

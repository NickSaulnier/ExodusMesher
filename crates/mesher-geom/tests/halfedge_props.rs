//! Half-edge connectivity against brute-force scans of the triangle list, on
//! small index ranges so non-manifold, inconsistent and degenerate edges are common.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_geom::{EdgeKind, HalfEdge, SurfaceMesh, TriangleSoup};
use proptest::prelude::*;

fn soup() -> impl Strategy<Value = (usize, Vec<[u32; 3]>)> {
    (1usize..8).prop_flat_map(|n| {
        (
            Just(n),
            prop::collection::vec(prop::array::uniform3(0..n as u32), 0..16),
        )
    })
}

fn ends(tris: &[[u32; 3]], h: u32) -> (u32, u32) {
    let t = tris[(h / 3) as usize];
    let c = (h % 3) as usize;
    (t[c], t[(c + 1) % 3])
}

fn brute_fans(tris: &[[u32; 3]], v: u32) -> usize {
    let inc: Vec<usize> = (0..tris.len()).filter(|&t| tris[t].contains(&v)).collect();
    let mut seen = vec![false; inc.len()];
    let mut fans = 0;
    for s in 0..inc.len() {
        if seen[s] {
            continue;
        }
        fans += 1;
        let mut stack = vec![s];
        seen[s] = true;
        while let Some(i) = stack.pop() {
            for j in 0..inc.len() {
                let shares_other = tris[inc[i]]
                    .iter()
                    .any(|&w| w != v && tris[inc[j]].contains(&w));
                if !seen[j] && shares_other {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
    }
    fans
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    #[test]
    fn matches_brute_force((n, tris) in soup()) {
        let pts = vec![[0.0; 3]; n];
        let mesh = SurfaceMesh::new(&TriangleSoup::new(&pts, &tris).unwrap());
        let hs = tris.len() as u32 * 3;

        // Edges: exactly the distinct unordered pairs, each with all its half-edges.
        let mut keys: Vec<[u32; 2]> = (0..hs)
            .map(|h| { let (a, b) = ends(&tris, h); [a.min(b), a.max(b)] })
            .collect();
        keys.sort();
        keys.dedup();
        prop_assert_eq!(mesh.edge_count() as usize, keys.len());
        for e in mesh.edges() {
            let [a, b] = mesh.edge_vertices(e).unwrap();
            let on_edge: Vec<u32> = (0..hs)
                .filter(|&h| { let (x, y) = ends(&tris, h); (x.min(y), x.max(y)) == (a, b) })
                .collect();
            prop_assert_eq!(mesh.edge_half_edges(e), &on_edge[..]);
            let kind = if a == b {
                EdgeKind::Degenerate
            } else if on_edge.len() == 1 {
                EdgeKind::Boundary
            } else if on_edge.len() == 2 {
                if ends(&tris, on_edge[0]).0 == ends(&tris, on_edge[1]).0 {
                    EdgeKind::Inconsistent
                } else {
                    EdgeKind::Manifold
                }
            } else {
                EdgeKind::NonManifold
            };
            prop_assert_eq!(mesh.edge_kind(e), Some(kind));
            for &h in &on_edge {
                let twin = (on_edge.len() == 2).then(|| HalfEdge(on_edge[usize::from(on_edge[0] == h)]));
                prop_assert_eq!(mesh.twin(HalfEdge(h)), twin);
                prop_assert_eq!(mesh.edge(HalfEdge(h)), Some(e));
            }
        }

        for v in 0..n as u32 {
            let out: Vec<u32> = (0..hs).filter(|&h| ends(&tris, h).0 == v).collect();
            prop_assert_eq!(mesh.outgoing(v), &out[..]);
            let incident: Vec<u32> = (0..tris.len() as u32).filter(|&t| tris[t as usize].contains(&v)).collect();
            prop_assert_eq!(mesh.vertex_triangles(v), incident);
            prop_assert_eq!(mesh.vertex_fan_count(v), brute_fans(&tris, v));
        }
    }
}

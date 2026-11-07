//! The volume mesh handed to the writers: points, tetrahedra grouped into
//! element blocks (one per part), and side sets.
//!
//! Conventions (shared by every writer):
//! - Tetrahedra are positively oriented: `orient3d(n0, n1, n2, n3)` is positive.
//! - Elements are numbered globally block by block, in block order: the
//!   first tet of block 2 follows the last tet of block 1.
//! - Sides use Exodus II TETRA4 side ordinals 1–4 ([`TET_SIDES`]).

use mesher_geom::Point3;
use mesher_geom::predicates::orient3d;
use mesher_geom::{Sign, vec3};

/// Exodus II TETRA4 sides: side `s` (1-based) has local nodes
/// `TET_SIDES[s - 1]` (0-based), listed so the face normal points out of a
/// positively oriented tetrahedron. In the 1-based numbering of the Exodus II
/// manual: side 1 = (1, 2, 4), side 2 = (2, 3, 4), side 3 = (1, 4, 3),
/// side 4 = (1, 3, 2).
pub const TET_SIDES: [[usize; 3]; 4] = [[0, 1, 3], [1, 2, 3], [0, 3, 2], [0, 2, 1]];

/// The local node opposite each side.
pub const TET_SIDE_OPPOSITE: [usize; 4] = [2, 0, 1, 3];

/// The nodes of side `side` (1-based) of `tet`, outward-ordered, or `None`
/// for a side number outside 1–4.
#[must_use]
pub fn tet_side_nodes(tet: [u32; 4], side: u8) -> Option<[u32; 3]> {
    let local = TET_SIDES.get(usize::from(side).checked_sub(1)?)?;
    Some(local.map(|i| tet[i]))
}

/// An element block: one part's tetrahedra.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
    pub name: String,
    pub tets: Vec<[u32; 4]>,
}

/// A side set: element sides, as `(global element, side)` with the element
/// 0-based in the global numbering and the side 1-based (Exodus ordinal).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SideSet {
    pub name: String,
    pub sides: Vec<(u32, u8)>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct VolumeMesh {
    pub points: Vec<Point3>,
    pub blocks: Vec<Block>,
    pub side_sets: Vec<SideSet>,
}

/// Something wrong with a mesh handed to a writer (a bug in the stage that
/// built it, not bad user input).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshProblem {
    NodeOutOfRange { element: u32, node: u32 },
    NotPositive { element: u32 },
    BadSide { set: String, element: u32, side: u8 },
    NonFinitePoint { node: u32 },
}

impl VolumeMesh {
    #[must_use]
    pub fn element_count(&self) -> usize {
        self.blocks.iter().map(|b| b.tets.len()).sum()
    }

    /// Every tet in global element order.
    pub fn elements(&self) -> impl Iterator<Item = [u32; 4]> + '_ {
        self.blocks.iter().flat_map(|b| b.tets.iter().copied())
    }

    /// Global number of the first element of each block.
    #[must_use]
    pub fn block_offsets(&self) -> Vec<u32> {
        let mut offsets = Vec::with_capacity(self.blocks.len());
        let mut next = 0u32;
        for b in &self.blocks {
            offsets.push(next);
            next = next.saturating_add(b.tets.len() as u32);
        }
        offsets
    }

    /// The faces of block `b` that no other tet of the block shares, as
    /// `(global element, side)`, in element order.
    #[must_use]
    pub fn boundary_sides(&self, b: usize) -> Vec<(u32, u8)> {
        let Some(block) = self.blocks.get(b) else {
            return Vec::new();
        };
        let offset = self.block_offsets().get(b).copied().unwrap_or(0);
        let mut faces: Vec<([u32; 3], u32, u8)> = Vec::with_capacity(block.tets.len() * 4);
        for (e, &tet) in block.tets.iter().enumerate() {
            for side in 1..=4u8 {
                if let Some(mut key) = tet_side_nodes(tet, side) {
                    key.sort_unstable();
                    faces.push((key, offset + e as u32, side));
                }
            }
        }
        faces.sort_unstable();
        let mut out = Vec::new();
        let mut i = 0;
        while i < faces.len() {
            let j = (i..faces.len())
                .find(|&j| faces[j].0 != faces[i].0)
                .unwrap_or(faces.len());
            if j - i == 1 {
                out.push((faces[i].1, faces[i].2));
            }
            i = j;
        }
        out.sort_unstable();
        out
    }

    /// Add a side set `<block>_boundary` for every block with a boundary.
    pub fn add_boundary_side_sets(&mut self) {
        for b in 0..self.blocks.len() {
            let sides = self.boundary_sides(b);
            if !sides.is_empty() {
                let name = format!("{}_boundary", self.blocks[b].name);
                self.side_sets.push(SideSet { name, sides });
            }
        }
    }

    /// Total tetrahedron volume.
    #[must_use]
    pub fn volume(&self) -> f64 {
        self.elements()
            .filter_map(|t| {
                let p = t.map(|v| self.points.get(v as usize).copied());
                let [Some(a), Some(b), Some(c), Some(d)] = p else {
                    return None;
                };
                let (u, v, w) = (vec3::sub(b, a), vec3::sub(c, a), vec3::sub(d, a));
                Some(vec3::dot(u, vec3::cross(v, w)) / 6.0)
            })
            .sum()
    }

    /// Check the conventions above. Returns the first problems found (at most
    /// `limit`).
    #[must_use]
    pub fn validate(&self, limit: usize) -> Vec<MeshProblem> {
        let mut problems = Vec::new();
        let n = self.points.len() as u32;
        for (i, p) in self.points.iter().enumerate() {
            if !p.iter().all(|c| c.is_finite()) {
                problems.push(MeshProblem::NonFinitePoint { node: i as u32 });
            }
        }
        for (e, tet) in self.elements().enumerate() {
            if let Some(&node) = tet.iter().find(|&&v| v >= n) {
                problems.push(MeshProblem::NodeOutOfRange {
                    element: e as u32,
                    node,
                });
                continue;
            }
            let [a, b, c, d] = tet.map(|v| self.points[v as usize]);
            if orient3d(a, b, c, d) != Sign::Positive {
                problems.push(MeshProblem::NotPositive { element: e as u32 });
            }
        }
        let count = self.element_count() as u32;
        for set in &self.side_sets {
            for &(element, side) in &set.sides {
                if element >= count || !(1..=4).contains(&side) {
                    problems.push(MeshProblem::BadSide {
                        set: set.name.clone(),
                        element,
                        side,
                    });
                }
            }
        }
        problems.truncate(limit);
        problems
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNIT: [Point3; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    /// The ordinals match the Exodus II manual's 1-based table exactly.
    #[test]
    fn side_ordinals_match_the_exodus_table() {
        let one_based: Vec<[usize; 3]> = TET_SIDES.iter().map(|s| s.map(|i| i + 1)).collect();
        assert_eq!(one_based, vec![[1, 2, 4], [2, 3, 4], [1, 4, 3], [1, 3, 2]]);
    }

    /// Every side, on a single positive tet, points away from the opposite node
    /// and covers the three other nodes.
    #[test]
    fn every_side_is_outward_on_a_single_tet() {
        assert_eq!(orient3d(UNIT[0], UNIT[1], UNIT[2], UNIT[3]), Sign::Positive);
        for side in 1..=4u8 {
            let nodes = tet_side_nodes([0, 1, 2, 3], side).unwrap();
            let opposite = TET_SIDE_OPPOSITE[usize::from(side) - 1] as u32;
            let mut all: Vec<u32> = nodes.to_vec();
            all.push(opposite);
            all.sort_unstable();
            assert_eq!(all, vec![0, 1, 2, 3], "side {side}");
            let [a, b, c] = nodes.map(|v| UNIT[v as usize]);
            // Outward: the opposite node is behind the face.
            assert_eq!(
                orient3d(a, b, c, UNIT[opposite as usize]),
                Sign::Negative,
                "side {side}"
            );
        }
        // Each side's face lies on the expected plane of the unit tet.
        let plane = |side: u8| tet_side_nodes([0, 1, 2, 3], side).unwrap();
        assert_eq!(plane(1), [0, 1, 3]); // y = 0
        assert_eq!(plane(3), [0, 3, 2]); // x = 0
        assert_eq!(plane(4), [0, 2, 1]); // z = 0
        assert_eq!(plane(2), [1, 2, 3]); // x + y + z = 1
        assert_eq!(tet_side_nodes([0, 1, 2, 3], 0), None);
        assert_eq!(tet_side_nodes([0, 1, 2, 3], 5), None);
    }

    #[test]
    fn boundary_of_two_tets_sharing_a_face() {
        let mut points = UNIT.to_vec();
        points.push([1.0, 1.0, 1.0]);
        // Second tet on the other side of face (1, 2, 3).
        let mesh = VolumeMesh {
            points,
            blocks: vec![Block {
                name: "b".into(),
                tets: vec![[0, 1, 2, 3], [1, 4, 2, 3]],
            }],
            side_sets: Vec::new(),
        };
        assert!(mesh.validate(10).is_empty(), "{:?}", mesh.validate(10));
        let sides = mesh.boundary_sides(0);
        assert_eq!(sides.len(), 6);
        assert!(!sides.contains(&(0, 2)));
        assert!((mesh.volume() - (1.0 / 6.0 + 1.0 / 3.0)).abs() < 1e-15);
    }

    #[test]
    fn validation_catches_inverted_and_bad_sides() {
        let mesh = VolumeMesh {
            points: UNIT.to_vec(),
            blocks: vec![Block {
                name: "b".into(),
                tets: vec![[0, 2, 1, 3], [0, 1, 2, 9]],
            }],
            side_sets: vec![SideSet {
                name: "s".into(),
                sides: vec![(0, 5)],
            }],
        };
        let problems = mesh.validate(10);
        assert!(problems.contains(&MeshProblem::NotPositive { element: 0 }));
        assert!(problems.contains(&MeshProblem::NodeOutOfRange {
            element: 1,
            node: 9
        }));
        assert_eq!(problems.len(), 3);
    }
}

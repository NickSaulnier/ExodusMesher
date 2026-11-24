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
    NodeOutOfRange {
        element: u32,
        node: u32,
    },
    NotPositive {
        element: u32,
    },
    BadSide {
        set: String,
        element: u32,
        side: u8,
    },
    NonFinitePoint {
        node: u32,
    },
    /// A face shared by more than two tets (one of them, and its side).
    FaceOnManyCells {
        element: u32,
        side: u8,
    },
}

/// A boundary patch of a [`FaceMesh`]: faces `start..start + count`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Patch {
    pub name: String,
    pub start: usize,
    pub count: usize,
}

/// A cell zone of a [`FaceMesh`] (an element block): cells
/// `start..start + count`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Zone {
    pub name: String,
    pub start: usize,
    pub count: usize,
}

/// A volume mesh as faces, in OpenFOAM's layout: the internal faces first,
/// ordered by owner and then neighbour (the owner the lower-numbered cell),
/// then each patch's faces together, patches in side set order. Each face is
/// listed out of its owner (its right-hand normal points from the owner
/// towards the neighbour, or out of the mesh); cells are the global
/// elements.
///
/// Every boundary face is in exactly one patch: the first side set listing
/// it other than a block's `<block>_boundary` set (an OBJ group's), else its
/// `<block>_boundary` set, which so holds the boundary faces in no group;
/// faces in no side set go to a patch `boundary`. Patches left empty are
/// omitted. Zones are the blocks.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FaceMesh {
    pub faces: Vec<[u32; 3]>,
    /// Per face.
    pub owner: Vec<u32>,
    /// Per internal face (the first `neighbour.len()` faces).
    pub neighbour: Vec<u32>,
    pub patches: Vec<Patch>,
    pub zones: Vec<Zone>,
    pub cells: usize,
}

impl FaceMesh {
    /// The number of internal faces.
    #[must_use]
    pub fn internal(&self) -> usize {
        self.neighbour.len()
    }
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

    /// The mesh as faces, for the face-based formats (OpenFOAM, Fluent):
    /// each unique face once, or the first problem found (a face on more than
    /// two cells, a node or side out of range).
    pub fn faces(&self) -> Result<FaceMesh, MeshProblem> {
        let n = self.points.len() as u32;
        let elements: Vec<[u32; 4]> = self.elements().collect();
        let mut keys: Vec<([u32; 3], u32, u8)> = Vec::with_capacity(4 * elements.len());
        for (e, &tet) in elements.iter().enumerate() {
            if let Some(&node) = tet.iter().find(|&&v| v >= n) {
                return Err(MeshProblem::NodeOutOfRange {
                    element: e as u32,
                    node,
                });
            }
            for side in 1..=4u8 {
                if let Some(mut key) = tet_side_nodes(tet, side) {
                    key.sort_unstable();
                    keys.push((key, e as u32, side));
                }
            }
        }
        keys.sort_unstable();

        // Each boundary face's patch: the first side set listing it other
        // than a block's `<block>_boundary`, else that set.
        let catch_all: Vec<bool> = self
            .side_sets
            .iter()
            .map(|s| {
                self.blocks
                    .iter()
                    .any(|b| s.name == format!("{}_boundary", b.name))
            })
            .collect();
        let mut listed: Vec<((u32, u8), usize)> = Vec::new();
        for pass in [false, true] {
            for (i, set) in self.side_sets.iter().enumerate() {
                if catch_all[i] != pass {
                    continue;
                }
                for &(element, side) in &set.sides {
                    if element as usize >= elements.len() || !(1..=4).contains(&side) {
                        return Err(MeshProblem::BadSide {
                            set: set.name.clone(),
                            element,
                            side,
                        });
                    }
                    listed.push(((element, side), i));
                }
            }
        }
        // The first listing wins: a stable sort keeps the pass order.
        listed.sort_by_key(|(k, _)| *k);
        listed.dedup_by_key(|(k, _)| *k);

        let mut internal: Vec<(u32, u32, [u32; 3])> = Vec::new();
        // Per side set, then for faces in none: the patch's faces.
        let mut by_set: Vec<Vec<(u32, [u32; 3])>> = vec![Vec::new(); self.side_sets.len() + 1];
        let mut i = 0;
        while i < keys.len() {
            let j = (i..keys.len())
                .find(|&j| keys[j].0 != keys[i].0)
                .unwrap_or(keys.len());
            let side_nodes = |k: usize| {
                let (_, e, s) = keys[k];
                (
                    e,
                    tet_side_nodes(elements[e as usize], s).unwrap_or(keys[k].0),
                )
            };
            match j - i {
                1 => {
                    let (e, nodes) = side_nodes(i);
                    let set = listed
                        .binary_search_by_key(&(e, keys[i].2), |(k, _)| *k)
                        .map_or(self.side_sets.len(), |at| listed[at].1);
                    by_set[set].push((e, nodes));
                }
                2 => {
                    // Owned by the lower cell, oriented out of it.
                    let (owner, nodes) = side_nodes(i);
                    internal.push((owner, keys[i + 1].1, nodes));
                }
                _ => {
                    return Err(MeshProblem::FaceOnManyCells {
                        element: keys[i].1,
                        side: keys[i].2,
                    });
                }
            }
            i = j;
        }
        internal.sort_unstable_by_key(|&(o, nb, _)| (o, nb));

        let mut fm = FaceMesh {
            faces: Vec::with_capacity(internal.len()),
            owner: Vec::with_capacity(internal.len()),
            neighbour: Vec::with_capacity(internal.len()),
            patches: Vec::new(),
            zones: Vec::new(),
            cells: elements.len(),
        };
        for (o, nb, nodes) in internal {
            fm.faces.push(nodes);
            fm.owner.push(o);
            fm.neighbour.push(nb);
        }
        for (set, mut faces) in by_set.into_iter().enumerate() {
            if faces.is_empty() {
                continue;
            }
            faces.sort_unstable_by_key(|&(e, nodes)| (e, nodes));
            let name = self
                .side_sets
                .get(set)
                .map_or_else(|| "boundary".to_owned(), |s| s.name.clone());
            fm.patches.push(Patch {
                name,
                start: fm.faces.len(),
                count: faces.len(),
            });
            for (e, nodes) in faces {
                fm.faces.push(nodes);
                fm.owner.push(e);
            }
        }
        let offsets = self.block_offsets();
        for (b, block) in self.blocks.iter().enumerate() {
            let start = offsets.get(b).copied().unwrap_or(0) as usize;
            fm.zones.push(Zone {
                name: block.name.clone(),
                start,
                count: block.tets.len(),
            });
        }
        Ok(fm)
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

    /// The right-hand normal of face `f` points away from node `away`.
    fn points_away(p: &[Point3], f: [u32; 3], away: Point3) -> bool {
        let [a, b, c] = f.map(|v| p[v as usize]);
        let n = vec3::cross(vec3::sub(b, a), vec3::sub(c, a));
        vec3::dot(n, vec3::sub(away, a)) < 0.0
    }

    fn centroid(p: &[Point3], t: [u32; 4]) -> Point3 {
        let q = t.map(|v| p[v as usize]);
        [0, 1, 2].map(|k| q.iter().map(|x| x[k]).sum::<f64>() / 4.0)
    }

    #[test]
    fn a_single_tet_is_four_boundary_faces_out_of_it() {
        let mut mesh = VolumeMesh {
            points: UNIT.to_vec(),
            blocks: vec![Block {
                name: "b".into(),
                tets: vec![[0, 1, 2, 3]],
            }],
            side_sets: Vec::new(),
        };
        mesh.add_boundary_side_sets();
        let f = mesh.faces().unwrap();
        assert_eq!((f.internal(), f.faces.len(), f.cells), (0, 4, 1));
        assert_eq!(
            f.patches,
            vec![Patch {
                name: "b_boundary".into(),
                start: 0,
                count: 4
            }]
        );
        let c = centroid(&mesh.points, [0, 1, 2, 3]);
        for face in &f.faces {
            assert!(points_away(&mesh.points, *face, c), "{face:?}");
        }
        assert_eq!(f.owner, vec![0; 4]);
    }

    /// Two tets sharing a face: one internal face, owned by the first and
    /// pointing into the second; an OBJ group's side set takes its faces out
    /// of `<block>_boundary`, which keeps the rest.
    #[test]
    fn two_tets_have_one_internal_face_and_patches_partition_the_boundary() {
        let mut points = UNIT.to_vec();
        points.push([1.0, 1.0, 1.0]);
        let tets = vec![[0, 1, 2, 3], [1, 4, 2, 3]];
        let mut mesh = VolumeMesh {
            points,
            blocks: vec![Block {
                name: "b".into(),
                tets: tets.clone(),
            }],
            side_sets: Vec::new(),
        };
        mesh.add_boundary_side_sets();
        // Tet 1's sides other than the shared one (its side 3, nodes 1, 3, 2).
        let group: Vec<(u32, u8)> = vec![(1, 1), (1, 2), (1, 4)];
        mesh.side_sets.push(SideSet {
            name: "top".into(),
            sides: group,
        });
        let f = mesh.faces().unwrap();
        assert_eq!(f.internal(), 1);
        assert_eq!((f.owner[0], f.neighbour[0]), (0, 1));
        let mut shared = f.faces[0];
        shared.sort_unstable();
        assert_eq!(shared, [1, 2, 3]);
        assert!(points_away(
            &mesh.points,
            f.faces[0],
            centroid(&mesh.points, tets[0])
        ));
        assert!(!points_away(
            &mesh.points,
            f.faces[0],
            centroid(&mesh.points, tets[1])
        ));
        // Patches in side set order, together covering the 6 boundary faces once.
        let names: Vec<(&str, usize)> = f
            .patches
            .iter()
            .map(|p| (p.name.as_str(), p.count))
            .collect();
        assert_eq!(names, vec![("b_boundary", 3), ("top", 3)]);
        assert_eq!(f.patches[0].start, 1);
        assert_eq!(f.patches[1].start, 4);
        for p in &f.patches {
            for k in p.start..p.start + p.count {
                let o = f.owner[k] as usize;
                assert_eq!(o, usize::from(p.name == "top"));
                assert!(points_away(
                    &mesh.points,
                    f.faces[k],
                    centroid(&mesh.points, tets[o])
                ));
            }
        }
        let mut all: Vec<[u32; 3]> = f
            .faces
            .iter()
            .map(|x| {
                let mut x = *x;
                x.sort_unstable();
                x
            })
            .collect();
        all.sort_unstable();
        all.dedup();
        assert_eq!(all.len(), 7);
        assert_eq!(
            f.zones,
            vec![Zone {
                name: "b".into(),
                start: 0,
                count: 2
            }]
        );
    }

    #[test]
    fn a_face_on_three_tets_is_a_problem() {
        let mut points = UNIT.to_vec();
        points.push([1.0, 1.0, 1.0]);
        points.push([2.0, 2.0, 2.0]);
        let mesh = VolumeMesh {
            points,
            blocks: vec![Block {
                name: "b".into(),
                tets: vec![[0, 1, 2, 3], [1, 4, 2, 3], [1, 5, 2, 3]],
            }],
            side_sets: Vec::new(),
        };
        assert!(matches!(
            mesh.faces(),
            Err(MeshProblem::FaceOnManyCells { .. })
        ));
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

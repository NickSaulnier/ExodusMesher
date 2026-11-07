//! Disjoint sets whose representative is always the smallest member, so
//! cluster and component ids are deterministic.

#[derive(Debug, Clone)]
pub(crate) struct DisjointSets {
    parent: Vec<u32>,
}

impl DisjointSets {
    pub(crate) fn new(n: usize) -> DisjointSets {
        DisjointSets {
            parent: (0..n as u32).collect(),
        }
    }

    pub(crate) fn find(&mut self, mut i: u32) -> u32 {
        while self.parent[i as usize] != i {
            let up = self.parent[self.parent[i as usize] as usize];
            self.parent[i as usize] = up;
            i = up;
        }
        i
    }

    pub(crate) fn union(&mut self, a: u32, b: u32) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra.max(rb) as usize] = ra.min(rb);
        }
    }
}

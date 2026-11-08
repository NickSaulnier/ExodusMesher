//! Opt-in timing for profiling the kernels. With the `profile` feature,
//! `prof::scope("name")` returns a guard that adds its lifetime to a
//! per-thread total for that name, and `prof::report()` lists the totals;
//! without it both compile to nothing. Timing never influences results.

#[cfg(feature = "profile")]
mod imp {
    use std::cell::RefCell;
    use std::collections::BTreeMap;
    use std::time::{Duration, Instant};

    thread_local! {
        static TOTALS: RefCell<BTreeMap<&'static str, (Duration, u64)>> =
            const { RefCell::new(BTreeMap::new()) };
    }

    pub struct Scope {
        name: &'static str,
        start: Instant,
    }

    impl Drop for Scope {
        fn drop(&mut self) {
            let elapsed = self.start.elapsed();
            TOTALS.with(|t| {
                let mut t = t.borrow_mut();
                let e = t.entry(self.name).or_default();
                e.0 += elapsed;
                e.1 += 1;
            });
        }
    }

    #[must_use]
    pub fn scope(name: &'static str) -> Scope {
        Scope {
            name,
            start: Instant::now(),
        }
    }

    /// The totals so far, slowest first, and clear them.
    #[must_use]
    pub fn report() -> String {
        TOTALS.with(|t| {
            let mut rows: Vec<(&'static str, (Duration, u64))> =
                std::mem::take(&mut *t.borrow_mut()).into_iter().collect();
            rows.sort_by_key(|r| std::cmp::Reverse(r.1.0));
            rows.iter()
                .map(|(name, (d, n))| {
                    format!(
                        "{:>10.3} s {:>10} calls {:>10.2} µs/call  {name}\n",
                        d.as_secs_f64(),
                        n,
                        d.as_secs_f64() * 1e6 / (*n).max(1) as f64
                    )
                })
                .collect()
        })
    }
}

#[cfg(not(feature = "profile"))]
mod imp {
    pub struct Scope;

    // Empty, so a scope can be ended early with `drop` in either build.
    impl Drop for Scope {
        #[inline(always)]
        fn drop(&mut self) {}
    }

    #[must_use]
    #[inline(always)]
    pub fn scope(_name: &'static str) -> Scope {
        Scope
    }

    #[must_use]
    pub fn report() -> String {
        String::new()
    }
}

pub use imp::{Scope, report, scope};

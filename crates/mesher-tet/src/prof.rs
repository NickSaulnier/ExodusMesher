//! Opt-in timing for profiling the kernels. With the `profile` feature,
//! `prof::scope("name")` returns a guard that adds its lifetime to a
//! per-thread total for that name, and `prof::report()` lists the totals;
//! without it both compile to nothing. Timing never influences results.
//! With `PROF_EVERY` set to a number of seconds, the totals so far also go
//! to stderr that often, for runs that may be stopped before they report.

#[cfg(feature = "profile")]
mod imp {
    use std::cell::{Cell, RefCell};
    use std::collections::BTreeMap;
    use std::sync::OnceLock;
    use std::time::{Duration, Instant};

    thread_local! {
        static TOTALS: RefCell<BTreeMap<&'static str, (Duration, u64)>> =
            const { RefCell::new(BTreeMap::new()) };
        static LAST_DUMP: Cell<Option<Instant>> = const { Cell::new(None) };
    }

    /// `PROF_EVERY`, read once.
    fn every() -> Option<Duration> {
        static EVERY: OnceLock<Option<Duration>> = OnceLock::new();
        *EVERY.get_or_init(|| {
            std::env::var("PROF_EVERY")
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|&s| s > 0.0)
                .map(Duration::from_secs_f64)
        })
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
                if let Some(every) = every() {
                    let now = Instant::now();
                    match LAST_DUMP.get() {
                        None => LAST_DUMP.set(Some(now)),
                        Some(last) if now - last >= every => {
                            LAST_DUMP.set(Some(now));
                            eprint!("--- profile so far\n{}", rows(&t));
                        }
                        Some(_) => {}
                    }
                }
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
        TOTALS.with(|t| rows(&std::mem::take(&mut *t.borrow_mut())))
    }

    fn rows(t: &BTreeMap<&'static str, (Duration, u64)>) -> String {
        let mut rows: Vec<(&'static str, (Duration, u64))> =
            t.iter().map(|(&k, &v)| (k, v)).collect();
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

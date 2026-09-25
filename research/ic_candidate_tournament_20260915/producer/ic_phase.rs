//! Exclusive instruction boundaries and native monotonic intervals.
//! Whole-process time is measured externally; the worker snapshot omits its
//! reporting tail. The evaluator charges that external remainder to setup.
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum Phase {
    Setup, FactorBase, Precompute, Queries, Pdp, RelationCheck, MatrixBuild,
    RelationLa, TargetDescent, RecoveryCheck, RhoSolve,
    TargetQuery, TargetPdp, TargetRelationCheck,
}

const PHASES: [Phase; 14] = [Phase::Setup, Phase::FactorBase, Phase::Precompute,
    Phase::Queries, Phase::Pdp, Phase::RelationCheck, Phase::MatrixBuild,
    Phase::RelationLa, Phase::TargetDescent, Phase::RecoveryCheck, Phase::RhoSolve,
    Phase::TargetQuery, Phase::TargetPdp, Phase::TargetRelationCheck];

impl Phase {
    fn label(self) -> &'static [u8] {
        match self {
            Self::Setup => b"ic_setup\0",
            Self::FactorBase => b"ic_factor_base\0",
            Self::Precompute => b"ic_precompute\0",
            Self::Queries => b"ic_queries\0",
            Self::Pdp => b"ic_pdp\0",
            Self::RelationCheck => b"ic_relation_check\0",
            Self::MatrixBuild => b"ic_matrix_build\0",
            Self::RelationLa => b"ic_relation_la\0",
            Self::TargetDescent => b"ic_target_descent\0",
            Self::RecoveryCheck => b"ic_recovery_check\0",
            Self::RhoSolve => b"reference_solve\0",
            Self::TargetQuery => b"ic_target_query\0",
            Self::TargetPdp => b"ic_target_pdp\0",
            Self::TargetRelationCheck => b"ic_target_relation_check\0",
        }
    }
    fn name(self) -> &'static str {
        let label = self.label();
        let name = std::str::from_utf8(&label[..label.len()-1]).expect("ASCII phase");
        name.strip_prefix("ic_").unwrap_or(name)
    }
}

struct Trace {
    phase: Phase,
    last: Option<Instant>,
    elapsed: [u64; 14],
    online_start: Option<Instant>,
    online_ns: Option<u64>,
}

impl Trace {
    const fn new() -> Self {
        Self { phase: Phase::Setup, last: None, elapsed: [0; 14],
            online_start: None, online_ns: None }
    }
}

thread_local! { static CURRENT: RefCell<Trace> = const { RefCell::new(Trace::new()) }; }

/// Enable native clocks at worker entry. Library unit tests need not enable them.
pub fn begin() {
    CURRENT.with(|slot| {
        let mut trace = Trace::new();
        trace.last = Some(Instant::now());
        *slot.borrow_mut() = trace;
    });
}

#[inline]
pub fn mark(next: Phase) {
    CURRENT.with(|slot| {
        let mut trace = slot.borrow_mut();
        if trace.phase != next {
            dump(trace.phase.label());
            if let Some(last) = trace.last {
                let boundary = Instant::now();
                let index = trace.phase as usize;
                trace.elapsed[index] += u64::try_from(boundary.duration_since(last).as_nanos()).expect("bounded interval");
                trace.last = Some(boundary);
            }
            trace.phase = next;
        }
    });
}

/// These endpoints exclude reusable IC preparation, fixture creation and report
/// serialization. Rho's own target-dependent setup belongs inside its interval.
pub fn begin_online(first: Phase) {
    mark(first);
    CURRENT.with(|slot| {
        let mut trace = slot.borrow_mut();
        assert!(trace.online_start.is_none() && trace.online_ns.is_none(), "one target per process");
        trace.online_start = Some(trace.last.expect("native clock enabled"));
    });
}

pub fn end_online() {
    mark(Phase::Setup);
    CURRENT.with(|slot| {
        let mut trace = slot.borrow_mut();
        let start = trace.online_start.take().expect("online interval started");
        trace.online_ns = Some(u64::try_from(trace.last.unwrap().duration_since(start).as_nanos()).expect("bounded interval"));
    });
}

/// Snapshot before constructing timing JSON. Its omitted tail and process
/// launch/exit remain charged by the external whole-process clock.
pub fn wall_snapshot() -> (BTreeMap<&'static str, u64>, Option<u64>) {
    CURRENT.with(|slot| {
        let trace = slot.borrow();
        let mut elapsed = trace.elapsed;
        if let Some(last) = trace.last {
            elapsed[trace.phase as usize] += u64::try_from(last.elapsed().as_nanos()).expect("bounded interval");
        }
        (PHASES.iter().map(|&p| (p.name(), elapsed[p as usize])).collect(), trace.online_ns)
    })
}

#[inline(never)]
fn dump(label: &'static [u8]) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        let args = [0x4354_0003usize, label.as_ptr() as usize, 0, 0, 0, 0];
        std::arch::asm!(
            "rol rdi, 3", "rol rdi, 13", "rol rdi, 61", "rol rdi, 51",
            "xchg rbx, rbx", in("rax") args.as_ptr(), inout("rdx") 0usize => _,
            out("rdi") _, options(nostack)
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = label;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn online_interval_closes_and_excludes_setup_on_both_sides() {
        begin();
        std::thread::sleep(Duration::from_millis(5));
        begin_online(Phase::TargetDescent);
        mark(Phase::TargetQuery);
        std::thread::sleep(Duration::from_millis(2));
        mark(Phase::TargetPdp);
        mark(Phase::TargetRelationCheck);
        mark(Phase::RecoveryCheck);
        end_online();
        std::thread::sleep(Duration::from_millis(5));
        let (wall, total) = wall_snapshot();
        let charged: u64 = ["target_descent", "target_query", "target_pdp", "target_relation_check", "recovery_check"].iter().map(|p| wall[p]).sum();
        assert_eq!(total, Some(charged));
        assert!(wall["setup"] >= 10_000_000);
        assert!(wall["target_query"] >= 2_000_000);
    }
}

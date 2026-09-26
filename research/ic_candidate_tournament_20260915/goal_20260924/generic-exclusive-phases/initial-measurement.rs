//! Opt-in exclusive phase clocks for the bounded generic IC worker.
//!
//! A session belongs to its controlling thread. A measured caller must keep
//! phase-changing work on that thread (a collector must use a serial outer loop)
//! and separately enforce its thread/CPU budget. Nested algebra may run under
//! the owning PDP interval. These labels deliberately differ from archived
//! producer schemas; a legacy interval parser must not silently admit them.
use serde::Serialize;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::rc::Rc;
use std::time::Instant;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(usize)]
pub enum Phase {
    Setup,
    FactorBase,
    Precompute,
    Queries,
    Pdp,
    RelationCheck,
    MatrixBuild,
    RelationLa,
    TargetQuery,
    TargetPdp,
    TargetRelationCheck,
    TargetDescent,
    RecoveryCheck,
    RhoSolve,
}

const PHASES: [Phase; 14] = [
    Phase::Setup,
    Phase::FactorBase,
    Phase::Precompute,
    Phase::Queries,
    Phase::Pdp,
    Phase::RelationCheck,
    Phase::MatrixBuild,
    Phase::RelationLa,
    Phase::TargetQuery,
    Phase::TargetPdp,
    Phase::TargetRelationCheck,
    Phase::TargetDescent,
    Phase::RecoveryCheck,
    Phase::RhoSolve,
];

impl Phase {
    fn label(self) -> &'static [u8] {
        match self {
            Self::Setup => b"generic_ic_setup\0",
            Self::FactorBase => b"generic_ic_factor_base\0",
            Self::Precompute => b"generic_ic_precompute\0",
            Self::Queries => b"generic_ic_queries\0",
            Self::Pdp => b"generic_ic_pdp\0",
            Self::RelationCheck => b"generic_ic_relation_check\0",
            Self::MatrixBuild => b"generic_ic_matrix_build\0",
            Self::RelationLa => b"generic_ic_relation_la\0",
            Self::TargetQuery => b"generic_ic_target_query\0",
            Self::TargetPdp => b"generic_ic_target_pdp\0",
            Self::TargetRelationCheck => b"generic_ic_target_relation_check\0",
            Self::TargetDescent => b"generic_ic_target_descent\0",
            Self::RecoveryCheck => b"generic_ic_recovery_check\0",
            Self::RhoSolve => b"generic_ic_rho_solve\0",
        }
    }

    fn name(self) -> &'static str {
        let label = self.label();
        std::str::from_utf8(&label[11..label.len() - 1]).expect("static ASCII phase")
    }

    fn online(self) -> bool {
        matches!(
            self,
            Self::TargetQuery
                | Self::TargetPdp
                | Self::TargetRelationCheck
                | Self::TargetDescent
                | Self::RecoveryCheck
                | Self::RhoSolve
        )
    }
}

struct Trace {
    generation: u64,
    phase: Phase,
    started: Instant,
    last: Instant,
    elapsed: [u64; 14],
    entered: [bool; 14],
    online_elapsed: [u64; 14],
    online_entered: [bool; 14],
    online_start: Option<Instant>,
    online_ns: Option<u64>,
    error: Option<&'static str>,
}

impl Trace {
    fn boundary(&mut self, next: Phase) {
        dump(self.phase.label());
        let now = Instant::now();
        let ns = u64::try_from(now.duration_since(self.last).as_nanos()).expect("bounded interval");
        self.elapsed[self.phase as usize] += ns;
        if self.online_start.is_some() {
            self.online_elapsed[self.phase as usize] += ns;
        }
        self.last = now;
        self.phase = next;
        self.entered[next as usize] = true;
        if self.online_start.is_some() && next.online() {
            self.online_entered[next as usize] = true;
        }
    }
}

#[derive(Default)]
struct State {
    generation: u64,
    trace: Option<Trace>,
}

thread_local! { static CURRENT: RefCell<State> = RefCell::new(State::default()); }

/// Raw observed intervals. The external runner must add launch, input before
/// session entry and the report/exit tail to setup; this is not a cold total.
#[derive(Debug, Serialize)]
pub struct Snapshot {
    pub schema_version: u8,
    /// An unentered phase stays unknown; the caller must explicitly justify
    /// absent/fused stages before an admission adapter can assign zero cost.
    pub phases_ns: BTreeMap<&'static str, Option<u64>>,
    pub observed_wall_ns: u64,
    pub online_phases_ns: BTreeMap<&'static str, Option<u64>>,
    pub online_wall_ns: Option<u64>,
}

/// A thread-bound guard. Errors and unwinding disable clocks automatically.
#[must_use]
pub struct Session {
    generation: u64,
    _owner: PhantomData<Rc<()>>,
}

impl Session {
    pub fn begin() -> Result<Self, &'static str> {
        CURRENT.with(|slot| {
            let mut state = slot.borrow_mut();
            if state.trace.is_some() {
                return Err("measurement session already active");
            }
            state.generation = state
                .generation
                .checked_add(1)
                .expect("bounded session count");
            let generation = state.generation;
            let now = Instant::now();
            let mut entered = [false; 14];
            entered[Phase::Setup as usize] = true;
            state.trace = Some(Trace {
                generation,
                phase: Phase::Setup,
                started: now,
                last: now,
                elapsed: [0; 14],
                entered,
                online_elapsed: [0; 14],
                online_entered: [false; 14],
                online_start: None,
                online_ns: None,
                error: None,
            });
            Ok(Self {
                generation,
                _owner: PhantomData,
            })
        })
    }

    pub fn finish(self) -> Result<Snapshot, &'static str> {
        CURRENT.with(|slot| {
            let mut state = slot.borrow_mut();
            let trace = state.trace.as_mut().ok_or("measurement session missing")?;
            if trace.generation != self.generation {
                return Err("measurement session changed");
            }
            if trace.online_start.is_some() {
                return Err("online interval not closed");
            }
            if let Some(error) = trace.error {
                return Err(error);
            }
            trace.boundary(Phase::Setup);
            let snapshot = Snapshot {
                schema_version: 1,
                phases_ns: PHASES
                    .iter()
                    .map(|&p| {
                        (
                            p.name(),
                            trace.entered[p as usize].then_some(trace.elapsed[p as usize]),
                        )
                    })
                    .collect(),
                observed_wall_ns: u64::try_from(
                    trace.last.duration_since(trace.started).as_nanos(),
                )
                .expect("bounded session"),
                online_phases_ns: PHASES
                    .iter()
                    .filter(|p| p.online())
                    .map(|&p| {
                        (
                            p.name(),
                            trace.online_entered[p as usize]
                                .then_some(trace.online_elapsed[p as usize]),
                        )
                    })
                    .collect(),
                online_wall_ns: trace.online_ns,
            };
            state.trace = None;
            Ok(snapshot)
        })
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        CURRENT.with(|slot| {
            let mut state = slot.borrow_mut();
            if state
                .trace
                .as_ref()
                .is_some_and(|t| t.generation == self.generation)
            {
                state.trace = None;
            }
        });
    }
}

#[inline]
pub fn enabled() -> bool {
    CURRENT.with(|slot| slot.borrow().trace.is_some())
}

/// Disabled callers neither read a clock nor issue a profiler request.
#[inline]
pub fn mark(next: Phase) {
    CURRENT.with(|slot| {
        let mut state = slot.borrow_mut();
        if let Some(trace) = state.trace.as_mut() {
            if trace.online_start.is_some() && !next.online() {
                trace.error = Some("non-target phase inside online interval");
            }
            if trace.phase != next {
                trace.boundary(next);
            }
        }
    });
}

/// One target per session. Preparation must be complete before this call.
pub fn begin_online(first: Phase) {
    CURRENT.with(|slot| {
        let mut state = slot.borrow_mut();
        if let Some(trace) = state.trace.as_mut() {
            if !first.online() || trace.online_start.is_some() || trace.online_ns.is_some() {
                trace.error = Some("invalid or repeated online start");
                return;
            }
            trace.boundary(first);
            trace.online_start = Some(trace.last);
            trace.online_entered[first as usize] = true;
        }
    });
}

pub fn end_online() {
    CURRENT.with(|slot| {
        let mut state = slot.borrow_mut();
        if let Some(trace) = state.trace.as_mut() {
            let Some(start) = trace.online_start else {
                trace.error = Some("online interval was not started");
                return;
            };
            trace.boundary(Phase::Setup);
            trace.online_start = None;
            trace.online_ns = Some(
                u64::try_from(trace.last.duration_since(start).as_nanos())
                    .expect("bounded online interval"),
            );
        }
    });
}

/// Return to the containing phase even on an early return from a substage.
#[must_use]
pub struct Scope {
    previous: Option<(u64, Phase)>,
    _owner: PhantomData<Rc<()>>,
}

pub fn scope(next: Phase) -> Scope {
    let previous = CURRENT.with(|slot| {
        slot.borrow()
            .trace
            .as_ref()
            .map(|t| (t.generation, t.phase))
    });
    mark(next);
    Scope {
        previous,
        _owner: PhantomData,
    }
}

/// Backend witness checking is shared by collection and individual descent.
pub fn relation_check_scope() -> Scope {
    let online = CURRENT.with(|slot| {
        slot.borrow()
            .trace
            .as_ref()
            .is_some_and(|t| t.online_start.is_some())
    });
    scope(if online {
        Phase::TargetRelationCheck
    } else {
        Phase::RelationCheck
    })
}

impl Drop for Scope {
    fn drop(&mut self) {
        if let Some((generation, previous)) = self.previous {
            let same = CURRENT.with(|slot| {
                slot.borrow()
                    .trace
                    .as_ref()
                    .is_some_and(|t| t.generation == generation)
            });
            if same {
                mark(previous);
            }
        }
    }
}

#[inline(never)]
fn dump(label: &'static [u8]) {
    #[cfg(test)]
    DUMPS.with(|count| count.set(count.get() + 1));
    #[cfg(target_arch = "x86_64")]
    unsafe {
        // Valgrind amd64 DUMP_STATS_AT: dump and reset, keep collection enabled.
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
thread_local! { static DUMPS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) }; }

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn exclusive_measurement_closes_online_and_observed_intervals() {
        let session = Session::begin().unwrap();
        mark(Phase::Precompute);
        std::thread::sleep(Duration::from_millis(1));
        begin_online(Phase::TargetQuery);
        std::thread::sleep(Duration::from_millis(1));
        {
            let _pdp = scope(Phase::TargetPdp);
            let _check = scope(Phase::TargetRelationCheck);
        }
        mark(Phase::RecoveryCheck);
        end_online();
        let snapshot = session.finish().unwrap();
        assert_eq!(
            snapshot.phases_ns.values().flatten().sum::<u64>(),
            snapshot.observed_wall_ns
        );
        assert_eq!(
            Some(snapshot.online_phases_ns.values().flatten().sum()),
            snapshot.online_wall_ns
        );
        assert!(snapshot.phases_ns["precompute"].unwrap() >= 1_000_000);
        assert!(snapshot.online_phases_ns["target_query"].unwrap() >= 1_000_000);
        assert_eq!(snapshot.phases_ns["relation_la"], None);
        assert_eq!(snapshot.online_phases_ns["rho_solve"], None);
        assert!(!enabled());
    }

    #[test]
    fn exclusive_measurement_disabled_is_noop_and_errors_drop_the_session() {
        let before = DUMPS.with(|count| count.get());
        mark(Phase::Pdp);
        begin_online(Phase::TargetQuery);
        end_online();
        drop(scope(Phase::RelationLa));
        drop(relation_check_scope());
        assert_eq!(before, DUMPS.with(|count| count.get()));
        let session = Session::begin().unwrap();
        assert!(Session::begin().is_err());
        begin_online(Phase::TargetQuery);
        assert_eq!(session.finish().unwrap_err(), "online interval not closed");
        assert!(!enabled());
        let session = Session::begin().unwrap();
        begin_online(Phase::TargetQuery);
        mark(Phase::Precompute);
        end_online();
        assert_eq!(
            session.finish().unwrap_err(),
            "non-target phase inside online interval"
        );
        assert!(!enabled());
    }

    #[test]
    fn exclusive_measurement_guard_is_thread_local_and_survives_unwind() {
        let session = Session::begin().unwrap();
        std::thread::spawn(|| {
            assert!(!enabled());
            mark(Phase::Pdp);
        })
        .join()
        .unwrap();
        let old_scope = scope(Phase::Precompute);
        drop(session);
        let next = Session::begin().unwrap();
        drop(old_scope); // Cannot change the next session's phase.
        assert_eq!(
            CURRENT.with(|s| s.borrow().trace.as_ref().unwrap().phase),
            Phase::Setup
        );
        drop(next);
        let _ = std::panic::catch_unwind(|| {
            let _session = Session::begin().unwrap();
            panic!("intentional unwind control");
        });
        assert!(!enabled());
    }

    #[test]
    fn exclusive_measurement_unstarted_and_repeated_online_are_distinct() {
        let snapshot = Session::begin().unwrap().finish().unwrap();
        assert_eq!(snapshot.online_wall_ns, None);
        assert!(snapshot.online_phases_ns.values().all(Option::is_none));
        assert_eq!(snapshot.phases_ns["pdp"], None);
        let session = Session::begin().unwrap();
        begin_online(Phase::TargetQuery);
        end_online();
        begin_online(Phase::TargetQuery);
        assert_eq!(
            session.finish().unwrap_err(),
            "invalid or repeated online start"
        );
        assert!(!enabled());
    }

    #[test]
    fn exclusive_measurement_relation_check_uses_the_correct_interval() {
        let session = Session::begin().unwrap();
        {
            let _check = relation_check_scope();
        }
        begin_online(Phase::TargetPdp);
        {
            let _check = relation_check_scope();
        }
        end_online();
        let snapshot = session.finish().unwrap();
        assert!(snapshot.phases_ns["relation_check"].is_some());
        assert!(snapshot.online_phases_ns["target_relation_check"].is_some());
        assert_eq!(
            Some(snapshot.online_phases_ns.values().flatten().sum()),
            snapshot.online_wall_ns
        );
    }
}

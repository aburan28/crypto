//! Exclusive Callgrind boundaries for the single-threaded optimized producer.
//! Startup, reporting and termination are setup. Descent does not enter the
//! relation-collection markers. No phase count is inferred from elapsed time.
use std::cell::Cell;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    Setup, FactorBase, Precompute, Queries, Pdp, RelationCheck, MatrixBuild,
    RelationLa, TargetDescent, RecoveryCheck, RhoSolve,
}

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
        }
    }
}

thread_local! {
    static CURRENT: Cell<Phase> = const { Cell::new(Phase::Setup) };
}

/// Finish the old interval before changing its attribution. The complete
/// process remains instrumented; marker overhead is charged, never subtracted.
#[inline]
pub fn mark(next: Phase) {
    CURRENT.with(|current| {
        let previous = current.get();
        if previous != next {
            dump(previous.label());
            current.set(next);
        }
    });
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

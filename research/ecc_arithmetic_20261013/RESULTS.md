# ECC arithmetic instruction-set pass: result

**Status:** measured field-level candidate; no qualified wall-time result.
This is a CPU binary-field implementation step within the larger requested
portfolio in [PLAN.md](PLAN.md). The exact run, controls and raw observations
are in the [run report](reports/run-20261013T1837Z/report.md), with
[architecture choices](ARCHITECTURE.md) and the
[open obstruction](../../problems/ecc-arithmetic-platform-measurement.md).

## Implementation

The ARM PMULL and Intel PCLMUL paths now use three carryless products for
two-limb field elements and six for three-limb field elements. The prior
schoolbook path used four and nine, respectively. The four-limb path and
portable fallback keep their earlier algorithm. Reduction, public types and
curve formulas are unchanged. This is an arithmetic operation-count change;
no claim is made that point or scalar arithmetic is faster. The portable
fallback still has operand-dependent control flow and table access, so this
change does not establish a constant-time secret-scalar path.

The modified public multiply passed 17 focused release tests, including
bitwise reference equivalence, boundary fields, reduction and inversion. An
x86-64 Linux cross-check of the frozen baseline and candidate compiled with
Rust 1.98.0; it did not execute on Intel hardware. The full release library
suite compiled and ran in the sandbox: 4,044 passed, nine unrelated local-socket
tests failed with `Operation not permitted`, and 114 were ignored. The
collaboration test group, including all nine, then passed 37/37 with loopback
access. A full-suite rerun with loopback access was stopped after hours under
host load above 150. There is no single completed full-suite pass.

## Measurement decision

The final source-bound local panel completed 24 paired AB/BA rounds per field
for GF(2^113), GF(2^131), GF(2^163) and GF(2^233), 500,000 public multiplies
per arm per round. All outputs matched the bitwise reference. GF(2^131)
baseline rounds ranged from 27,826,500 to 932,235,083 ns; candidate rounds
ranged from 16,801,958 to 694,301,375 ns. The unchanged four-limb control
also varied substantially. The host load was above 100 during the panel.
Those observations fail the controlled wall-time comparison gate. There is
no admitted speedup or regression from these timings.

![Every final paired round](paired_mac_final.svg)

No Intel or GPU job completed. A sandboxed Modal attempt could not connect;
automatic approval review rejected the network-enabled retry because it
would export the source bundle to Modal without specific authorization.
The user can approve that exact source export after reviewing the bundle, or
an already-authorized worker can be used. The broader prime curves,
Curve25519/Ed25519, Koblitz scalar path, Metal and CUDA/HIP kernels remain
unimplemented or unmeasured as listed in [PLAN.md](PLAN.md).

The ECC2K boundary scoreboard and ecbench comparisons were checked at base
commit `11459b07f` and receive no updated performance claim, because only
a noisy field subroutine panel completed. The study-local graph is the visual
record of this measured scope.

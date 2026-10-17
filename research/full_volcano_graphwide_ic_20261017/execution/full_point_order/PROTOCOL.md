# Full selected-point-order C diagnostic

Frozen before measurement. This extends the separate order-7 diagnostic; it
does not rewrite it or replace the original full graph/method-matrix request.
Original C has p=1637958981148213249 and N=p-1=2^26*3^20*7. Both arms receive
the same public (curve,G,Q,M,factors), with certified exact ord(G)=M dividing N.
No planted scalar is an argument to either solver.

Use the immutable 1,869-vertex C prefix with SHA-256
468ce8320b6b7deae37198ca3afc07534fa4fe0e950f621e45ad322f23a2c881.
Require its native map replay receipt. Join existing local labels by exact
DB hash and vertex id: 517 exact vertices / 24 observed strata, 1,352 unknown.
Full C remains 1,934,901,248 vertices / 165 strata and is not completed here.

For each model, deterministically search at most 4,096 x coordinates. Compute
each candidate's exact point order by stripping the complete known prime
factorization of N and checking annihilation and every remaining prime divisor.
Choose the first point whose order contains 2,3,7. Keep its entire exact order;
do not cofactor-project it to a smaller benchmark subgroup. This certifies M,
not maximality or cyclicity of E(Fp). Generate two separate public known-answer
targets from a fixed hash seed modulo M. Keep scalar truth outside solver
inputs. Generator/order/target construction is shared fixture preparation,
recorded separately and excluded from both solver intervals.

Arm ec.ph_bsgs.full_point_v1 performs standard Pohlig--Hellman on all prime
powers in M, using a reusable bounded baby-step/giant-step digit table for
each prime, then CRT. This is an explicit direct control, not a fastest
external implementation or a generic rho run on the unfactored order.

Arm ff.k1_ph_relations.full_point_v1 projects to the same prime-power
components. Standard reduced Tate/Miller evaluates t_m(G_m,R) and
t_m(Q_m,R), varying the first argument. Search at most 256 auxiliary R points
per component and at most 128 shifted divisors per evaluation. Require the
generator image to have exact order m; reject degenerate transfers. No fake
pairing, discrete-log fallback or assumed injectivity is allowed.

Use multiplicative PH digit lifting in the destination. For primes 2,3,7,
build up to two distinct inverse classes of small-prime images; prime 2/3
use one class. Recover their initially unknown logs by checked one-/two-factor
relations and modular Gaussian elimination. The identity is a decomposition
slot, not a factor-base unknown; one-factor relations handle characteristic-2
digit coefficients correctly. Search at most q anchor and q target shifts.
Reuse verified base logs across digits of one component within one target.
Log collection queries, pair tests, equations, rank and matrix-entry work
separately for setup and target work. Verify base logs and every digit.

Both arms reconstruct modulo M by CRT. A runner check verifies [x]G=Q and the
separate public truth. Reject outside-subgroup aliases even when [M]Q=O.
Never suppress failed rows or substitute a successful control result.

Before execution, test small p=109/113/127 controls with exact enumeration,
all scalar inputs of selected control subgroups, both first-argument pairing
linearity and full image order, degenerate/invalid inputs, and original C
known-answer targets. Preserve all test failures. Then run source-only and
all-prefix populations separately, two one-target workloads per vertex and
both arms: 7,476 method records in the prefix run. Re-execute into fresh output
and compare every non-clock field natively. This is same-implementation replay,
not independent arithmetic proof.

Cold solver intervals include exact input-order certificate verification,
projection, all method setup, auxiliary-point
search, pairing, base/index/matrix work, digit recovery, CRT and final EC
verification. Input/output and shared fixture generation remain outside.
Record phase operation ledgers and unqualified end-to-end observations;
square-root internals, matrix inversions, allocation and integer work remain
explicitly unpriced. No timing speedup, fastest-method, scaling, sub-rho or
security claim is admitted. This is projected-subgroup relation algebra on
smooth known-answer fixtures, not NFS/integer-smoothness IC or a large-k
finite-field backend. The original eight-target/five-round matrix, transported
targets, independent replay and other original classes remain open.

Run through the repository isolation wrapper on CPU 4 using its unchanged
default admission thresholds. Preserve any refusal. If admission cannot be
earned, do not relax thresholds: count-only execution under the busy lock is
permitted with all wall observations explicitly unqualified. Finish all builds
and tests before execution; do not run other heavy work during measurement.

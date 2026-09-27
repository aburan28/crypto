# Frozen protocol: Frobenius-rotated normal-basis support gate

Status: source/input freeze after preregistration; no support or lift-density outcomes have been run.
This gate follows the complete toy four-sum oracle in
[PR #757](https://github.com/aburan28/crypto/pull/757). It is a bounded
representation and exact group-law support experiment, not
an S6/S7 solver, relation collection, rank, rho, or ECC2K-130 log claim.
It follows Galbraith–Granger–Merz–Petit, *On Index Calculus Algorithms for
Subfield Curves*, §3.1 (https://eprint.iacr.org/2020/1315.pdf). Their
normal-basis construction uses distinct x-spaces
V_i = span_{F2}{beta^(2^(m*j+i)): 0 <= j < d}, 0 <= i < m, m*d <= n,
and F_i = tau^i(F_0). Its m! relation-count advantage is conditional on
point-decomposition solve costs and cannot be read as an attack speedup.

## Frozen models and construction

Both models use E: y^2 + x*y = x^3 + 1 over F_(2^n), with full point
coordinates and the point at infinity. The exact n=131 structural/density
model uses polynomial x^131+x^13+x^2+x+1 (integer
0x800000000000000000000000000002007), subgroup order
q=680564733841876926932320129493409985129 and cofactor 4. This
polynomial-basis model is field-isomorphic to the public challenge field;
no Certicom target coordinate conversion or challenge log is attempted.
The exact small model uses n=13, polynomial x^13+x^4+x^3+x+1
(0x201b), and must independently count #E=8012=4*2003 with 2003 prime.
Before any support outcome, the implementation must verify both field
polynomials by the prime-degree Rabin irreducibility test. Both producer
and independent replay derive the group orders from #E(F_2)=4 and the
Weil-trace recurrence t_0=2, t_1=-1,
t_n=-t_(n-1)-2t_(n-2), #E(F_(2^n))=2^n+1-t_n; it returns
8012 for n=13 and 2722258935367507707729280517973639940516=4q
for n=131.

For each n, beta_A=3=1+x is accepted only if its n Frobenius conjugates
have exact F2 rank n and tau^n(beta_A)=beta_A. At n=13, beta_B is the
smallest integer >3 that is normal and not a Frobenius conjugate of
beta_A; the field-only preflight gives beta_B=7. It is chosen using only
field arithmetic, before examining any factor-base or target support. Both beta choices are measured, never
selected post hoc. For n=131 the fixed beta is beta_A.

At n=131 the cells are (m,d)=(5,24),(5,25),(5,26),(6,20),(6,21).
For each, verify rank(V_i)=d, V_i∩V_j={0} for i!=j, rank(sum V_i)=m*d,
and tau(V_i)=V_(i+1) for 0<=i<m-1 by exact F2 elimination.
The normal-basis trace must satisfy Tr(beta)=1 and
Tr(sum_j c_j beta^(2^(m*j+i))) = parity(c) for every tested mask.
For the n=131 compressed-column structural check, deterministically select
the first x>=2 with a rational lift whose [4]-projection H is nonzero;
construct the lift by half-trace. Assert [q]H=O, lambda^2+lambda+2=0
mod q, lambda^131=1 mod q and tau(H)=[lambda]H by actual group law for
lambda=196511074115861092422032515080945363956. Replay the selected
point, [4]tau(P)=tau[4](P), and scalar identity through separate field
and inverse arithmetic. The deterministic structural preflight selects
x=3, P=(3,780119811075450012506287188407015201238) and
H=[4]P=(302709522240084841455960327861177137188,
1704224067622899467506648521202877826049). This certifies the
action on the cyclic prime subgroup under the checked #E=4q model; no
target log is calculated. The lambda constant is consistent with the
direct public challenge-point identity in
`RESEARCH_ECC2K130_RELATION_SWEEPS.md` §4.3; this gate separately
checks it in the chosen polynomial basis.
The n=13 exact cells are m=5 or 6 with d=2 and the same construction.
Each F_i contains *all* rational point lifts of x in V_i, including
both signs when x!=0, the one point (0,1), and no infinity. The
equal-cardinality control repeats F_0 in every summand slot. No
post-result factor-base adjustment is allowed.

## Exact cofactor and support checks

Verify O,(0,1),(1,0),(1,1) as precisely the rational [4]-kernel:
(0,1) has order 2 and each x=1 point doubles to (0,1).
For every n=13 target Q=[k]H, 0<=k<2003, use the same subgroup
generator H=[4]G, where G is the first affine point in ascending
(x,y) order for which [4]G has order 2003. Verify tau(H)=[89]H by full
group law; 89 is the unique root of lambda^2+lambda+2=0 mod 2003
whose 13th power is 1. Verify [4]tau(P)=tau[4](P) for every point in
every constructed F_i, and map F_i to F_0 by tau powers.

For both rotated and repeated controls, enumerate all labelled
m-tuples by independent exact group-law algorithms. Record the full
point-sum histogram, exact distinct support, tuple collisions and
one witness for every supported sum. For each Q and each of the four
torsion points T, retain the exact multiplicity and witness of Q+T.
Also compute the [4]-projected multiset of sums and verify, for every Q,
that [4]S=[4]Q exactly when some S=Q+T. Record |F_i|, sign-pair
count, |[4]F_i|, projected duplicates, physical point choices
sum_i|F_i|, and compressed log columns from one projected F_0.
The latter compression is justified only after the lambda action
check; it is not a measured matrix rank. Report support over both all
2003 subgroup targets and 2002 nonzero targets. Independently replay
every positive witness and every negative support decision by a second
full tuple enumeration with separately implemented field inversion.

For every n=13 H target, retain machine-readable k, coordinates,
Tr(x_Q) (O separately), tau-orbit length, cofactor class, the four
Q+T multiplicities, projected multiplicity and support flags. Report
strata for those cheap public features, including features constant
on all nonzero H targets. No classifier is fit to this complete census;
any later classifier requires a new preregistered held-out gate.

## n=131 density and count admission

No n=131 factor base is enumerated. A separate covariance panel uses
256 SHA-256-derived coefficient masks per (m,d), with 0, 1, the top
single-bit and all-ones masks included as edge controls. Identical masks
are transported across all i by tau^i; each resulting x must have the
same curve-liftability as its source, and the trace parity must agree.
A separate density panel uses 2^14 distinct nonzero SHA-256-derived masks
per (m,d), excluding covariance masks. The domain is
ECC2K130-ROTATED-SUBSPACE-20260925-v1/{cov|density}/m/d/counter;
digest bytes are read big-endian and reduced modulo 2^d; duplicates,
zero and excluded masks are skipped by advancing the counter. Freeze
the generated mask-file hashes and beta before running density. Count
the trace-solvable x, assigning 2 rational lifts to each solvable
nonzero x and 0 otherwise; x=0 contributes exactly one. Check the
covariance panel independently across all i. Give Wilson intervals
only as descriptive pseudorandom-mask model summaries, not rigorous
bounds for the complete V_i. The density sample does not measure
n=131 target support or solve any PDP.

The pre-outcome count grid reports, under the explicit idealization
|F_i|≈2^d, the necessary raw full-curve tuple ceiling
min(1,2^(m*d)/(4q))≈min(1,2^(m*d-131)), and the projected
prime-subgroup ceiling min(1,2^(m*d)/q), at most four times the raw
ratio. Use exact q in the emitted table and distinguish this
conditional estimate from a bound using measured |F_i| or
|[4]F_i|. The 1% projected necessary threshold first clears at
m=5,d=25 and m=6,d=21; it promises no actual coverage.

| m | d | md | Idealized raw ceiling / (4q) | Idealized projected ceiling / q | Physical choices m·2^d | Compressed F0 proxy |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 24 | 120 | 0.00048828125 | 0.001953125 | 83,886,080 | 16,777,216 |
| 5 | 25 | 125 | 0.015625 | 0.0625 | 167,772,160 | 33,554,432 |
| 5 | 26 | 130 | 0.5 | 1 (ratio exceeds 1) | 335,544,320 | 67,108,864 |
| 6 | 20 | 120 | 0.00048828125 | 0.001953125 | 6,291,456 | 1,048,576 |
| 6 | 21 | 126 | 0.03125 | 0.125 | 12,582,912 | 2,097,152 |

The decimals display the exact-q formulas rounded to the shown precision.
The physical and compressed counts are idealized point-size proxies, not
measured matrix columns or upper bounds; rational lift density, [4] duplicates and signs
will be reported separately. The d=25/m=5 and d=21/m=6 materialized
three-summand projections at 24 bytes/tuple are approximately 906.7 ZB
and 221.4 EB respectively, conditional on no deduplication.

Physical point choices are roughly m*2^d, while Frobenius/known-lambda
compression can represent log columns through F_0. The displayed 2^d
compressed count is an idealized proxy, not an upper bound: with both
rational signs, |[4]F_0| <= |F_0| <= 2^(d+1)-1. Exact n=13 projected
duplicates and the n=131 point check must be reported before interpreting
this proxy as useful columns. A materialized
three-summand MITM half at d=25 or 21 conditionally holds 2^75 or
2^63 tuples before collisions, respectively; price 24 bytes/record
as a conditional storage projection, never a lower bound for a
different implicit solver.

## Decision, costs and preservation

Continue toward a new implicit S6/S7 solver-design gate only if n=13
rotated projected support is at least 1.25 times its same-beta repeated
control for one m, every exact check passes, and n=131 density/structure
does not falsify the count admission. If both m values have projected
support ratio <=1.10 on both beta choices, mark only the n13
*support-improvement* hypothesis negative; a separately budgeted
solver-symmetry gate remains possible if the exact n131 count admission
survives. Intermediate outcomes are inconclusive. Any (m,d) whose
sample-calibrated *necessary* projected tuple ceiling remains below
1% even under a descriptive 95% density upper interval is rejected
for the 1% admission target. A positive toy result only authorizes a
separate solver architecture/memory experiment; it does not establish
n=131 support, a relation rank, or a speed advantage.

Charge field multiplication, squaring, inversion, point addition and
scalar multiplication for field/normal-basis construction, full
point-lift setup, every tuple or dynamic-programming state, all four
Q+T right-hand sides per target even when one succeeds, [4] projection,
and full verification. Tuple multiplicity and distinct target support are
separate metrics throughout.
Record CPU, wall and sampled peak RSS as secondary. Bound each
small-model variant by 180 s and 512 MiB and each n=131 density cell
by 120 s and 512 MiB. On a cap or validation failure, preserve its
receipt and stop the affected cell; no result substitution or silent
target reduction. Archive source/input hashes, raw per-target counts,
all failed attempts, independent replay and the resulting decision in
one focused PR with CI and the canonical scoreboard/ledger update.

## Source/input freeze and exact runner

This section fixes the implementation and inputs before the first support or
lift-density measurement. The field-only Rabin/rank/trace preflight checked
n=13 polynomial `0x201b`, n=131 polynomial
`0x800000000000000000000000000002007`, beta_A=3 and first alternate
n=13 beta_B=7; it examined no support or density outcomes. The complete
machine-readable hash map is `FROZEN.json` (SHA-256
`77c313cc3db97430e77a797492c32037eecb6189225fc478789d8db6f659edc3`).
The input manifest SHA-256 is
`b1bfffd8ac22d47cd10f73e3745e0b295d6466b0c300df6674e7dd8d46381da4`.
The source SHA-256 values are:

| File | SHA-256 |
| --- | --- |
| `gate.py` | `d3f0f6e1515282a25efb4eb6a5f68ad754f1d9be68e4a741f17a748d9a58d782` |
| `verify.py` | `8b665a8d5a1d92cf17f64106dab4bf7134ecaffd06e669404a7508246282c1c8` |
| `run.py` | `052553f2dedacdaf23da6230bdfd1625c02bb3375718b6c67c087d086fa62618` |
| `ci_replay.py` | `f50bc43fc70b9976e04283c15d6fa2b93a84358bd1985374c2044fb47d455683` |

Use Python >=3.12 and run from the repository root:

```sh
python3 research/notes/ecc2k130/rotated_subspace_support_20260925/ci_replay.py
python3 research/notes/ecc2k130/rotated_subspace_support_20260925/run.py --out /private/tmp/rotated-subspace-run-20260925
```

`run.py` verifies the hash freeze before invoking the n13 toy producer,
n131 density producer and separate verifier, preserving stdout, stderr,
UTC timestamps, source/input hashes, command status and every raw file's
hash even on a run failure. The n13 producer charges one shared cold
curve/field/group/target setup plus stage-separated per-variant normal/factor
construction, exact histogram, `[4]` projection/witness self-check and all
four Q+T right-hand sides. The independent verifier has a separate charged
receipt. A counterfactual cold single-variant cost adds the shared setup to
that variant's stage total; the full eight-arm campaign charges it once.
The n131 density producer likewise separates shared field/normal/trace-mask
setup from each cell. `ru_maxrss` is reported as process high-water RSS,
converted to bytes; per-arm values can include a prior arm's high-water mark.
A POSIX SIGALRM wall deadline of 180 seconds for each toy arm and 120
seconds for each density cell covers its construction, enumeration and file
emission; an arm that expires is rejected and leaves a quantitative producer
failure file alongside already completed cells. The 512-MiB RSS cap is a
retrospective admission check after each cell, because portable Python
does not enforce a hard per-cell memory limit. An over-cap cell is rejected
and retained as a failed attempt. The aggregate runner timeouts are backup
process guards and do not replace these per-cell deadlines. Python 3.12 and 3.13 use the same integer arithmetic;
wall and CPU are host-specific diagnostics, while operation counts and exact
support are primary. The focused CI checks this freeze even before evidence
exists; after evidence is added it replays the entire archived census and
density sample with independent arithmetic.

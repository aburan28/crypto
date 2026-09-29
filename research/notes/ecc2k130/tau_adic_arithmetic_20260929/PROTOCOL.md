# Tau-adic known-scalar arithmetic diagnostic

Frozen 2026-09-29 UTC, before measured runs. This is a standalone arithmetic
study using synthetic points and known scalars. It contains no rho solver,
challenge target, walk replacement, collector, or cloud job.

## Hypothesis and boundaries

Reduced tau-NAF can reduce the cost of *known-scalar multiplication* on
`E_a: y^2 + xy = x^3 + a*x^2 + 1`, since tau replaces doubling with
coordinate squaring. Unreduced tau-NAF can be almost twice as long and must
be measured as a control. This does not establish a walk or ECDLP gain.

The reference is matched binary double-and-add, with binary NAF as a stronger
control, using the same affine arithmetic backend. Sage's native point
multiplication is also timed and independently verifies every answer.
Native Sage is a separate backend diagnostic, not the shipping GPU baseline.

For the whole application, if an eligible operation occupies fraction `f`
and its local cost ratio is `r = candidate/reference`, Amdahl's bound is
`total ratio = (1-f) + f*r`. The eligible fraction is unmeasured. Set GPU
throughput, whole-walk cost, ECDLP operation count, and end-to-end speedup to
null. Do not convert CPU arithmetic timings into GPU rates.

## Frozen inputs and accounting

- Audit repository commit `524d38773a5ee5ba8f7ffe7e969241f2218dc1e4`.
- Benchmark E_0 at m=31,83,131. Polynomial taps: (3,0), (45,2,1,0),
  (13,2,1,0), with the leading term z^m implicit. Check irreducibility.
- Training seed 20260929, independent holdout seed 20260930; 24 cases per
  panel and degree. Generate distinct synthetic points, clear cofactor 4,
  and choose known nonzero scalars below the order divided by 4.
- Exhaustively verify all E_0 and E_1 points over GF(2^5) with scalars
  -17 through 17, plus edge scalars 0,1,2,3,2^(m-1),2^m-1,2^m,2^m+1
  and their negatives on each large curve. Include infinity and torsion.
- Verify tau^2 - mu*tau + [2] = 0 for mu=2*a-1 and tau^m(P)=P.
- Reduce in Z[tau] modulo tau^m-1, which annihilates every rational point.
  Use exact integer arithmetic and a fixed 3x3 neighborhood of rounded
  quotient coefficients; no claim of globally minimal representatives.
- Compare binary, binary NAF, unreduced tau-NAF, reduced tau-NAF, and Sage
  native multiplication. Charge per-scalar recoding and reduction. Field/
  curve construction is shared fixture setup, timed separately and excluded
  from the per-scalar diagnostic. No precomputed digit tables.
- Warm each method. Record five A/A pairs for the binary-NAF control and
  seven A/B rounds with the method order reversed in alternate rounds.
  Record medians, minima, every sample, output hashes and operation vectors.
  Verification, hashing, and counters are outside the timed regions;
  counted reruns must match the same verified outputs.
- Pin one logical CPU; attempt matching NUMA memory placement and record
  the result. Pinning does not reserve a core. Record CPU family/model,
  virtualisation, visible NUMA topology, memory limits/type when available,
  Sage/Python versions, process/load snapshots, source/input/output hashes.
  Shared-VM timings are stage diagnostics with a contention limitation.

## Success and stop rules

Any correctness mismatch stops the experiment. A promising arithmetic result
requires reduced tau-NAF median cost at least 10% below binary NAF on both
m=83 and m=131 and on the independent holdout, with all answers equal. The
improvement must exceed the largest observed A/A deviation from parity.
These are screening criteria, not a statistically established GPU gain.

Do not change production kernels based on this study. If the audited hot
loop contains no long known-scalar multiplication, reject direct insertion
of tau-NAF into that loop as unsupported. An existing sparse Frobenius/add
expression is not a binary double-and-add operation awaiting recoding.

Stop after the frozen panel and correctness tests. Preserve failures and
regressions. Missing CUDA hardware/toolchain is a recorded limitation;
do not launch a distributed discrete-log computation.

## Reproduce

From this directory with SageMath 10.9:

```bash
sage -python benchmark.py --output results/run-01.json
```

The output path must not already exist. A failing run writes its exception
and partial evidence to that path before returning a nonzero status.

## References

- Ahmadi, Hankerson, Rodriguez-Henriquez, *Parallel Formulations of Scalar
  Multiplication on Koblitz Curves* (2007), sections 1-2:
  https://cacr.uwaterloo.ca/techreports/2007/cacr2007-18.pdf
- Avanzi, Heuberger, Prodinger, *Minimality of the Hamming Weight of the
  tau-NAF for Koblitz Curves*:
  https://www.math.aau.at/heuberger/publications/pdf/tauextcrypt.pdf

The exact quotient rounding used here is a simple reference construction,
not an implementation or performance claim for Solinas's partial reduction.

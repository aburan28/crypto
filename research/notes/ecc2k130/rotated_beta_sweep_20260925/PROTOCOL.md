# Preregistered n19 rotated normal-generator support sweep

Status: source and selection rule preregistered before any new six-summand
support outcome. This is a small-curve factor-base sensitivity experiment,
dependent on the merged exact point corpus in
[PR #767](https://github.com/aburan28/crypto/pull/767) (merge
`c77767c4a653734f428e110cca29985d721476b2`). Its archived
`evidence/raw.tar.gz` SHA-256 is
`39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`.
The reference is that archive's n19, beta-3, m6,d2 row: six 7-point factors,
117649 labelled tuples, 62389 projected supported subgroup targets among
`q=130873`, and 68484 misses. The source curve is
`E:y^2+xy=x^3+1` over `F_2^19` with polynomial `0x80027`,
`#E=523492=4q`, `H=(385982,301867)`, and Frobenius scalar `41811`.
No curve, polynomial, rung, arity, subgroup or parent corpus may be replaced
after the sweep observes support.

## Hypothesis, boundary and decision

At fixed n19 model, arity and factor size, distinct normal generators may
change projected six-sum support and which *fixed* #767 point targets have
representations. This is a sensitivity hypothesis, not an n131 prediction.
The immutable reference is beta 3 from #767; the necessary counting bound
for every admitted 7-point candidate is at least
`q - 7^6 = 13224` unsupported subgroup targets. The predeclared descriptive
effect threshold is an absolute support change of at least 1% of q (at least
1309 subgroup targets) relative to beta 3. A second, separately reported
discriminator is any exact membership flip among #767's eight n19 Q targets.
Both can fail; all zero/negative outcomes and failed arms stay in the PR.
Neither threshold is a speed boundary or an attack improvement.

The **follow-up priority gate**, distinct from the 1%-of-q descriptive
threshold, requires one admitted candidate at identical seven physical F0
points and three nonzero negation-paired projected columns (seven signed
projected points including O) to have at least 10% more projected support
than beta 3: `S_candidate >= ceil(1.10 * 62389) = 68628`. It must also have
no worse projected tuple energy *per supported target*:
`E_candidate / S_candidate <= E_beta3 / 62389`, checked by exact integer
cross multiplication, where `E = sum_R multiplicity(R)^2`. If no candidate
passes both conjuncts, prioritize other factor-base/solver work over further
beta tuning. This is a toy-rung decision rule, not evidence of transfer to
n131 or ECC2K-130.

The admission gate is at least three and at most four *preselected* additional
normal generators with complete producer and independent-verifier success
under the caps below. If the frozen selection rule yields fewer than three,
run exactly those admitted by its deterministic fallback and report the
admission failure; do not choose another beta after seeing any support result.
No solver, SAT/UNSAT engine, relation rank, logarithm, full-cost `S`, matched
rho ratio, or ECC2K-130 extrapolation is measured here.

## Selection rule frozen before support

The exact domain string is `ECC2K130-ROTATED-BETA-SWEEP-20260925-v1`.
For counters `0..4095`, compute `beta = 1 +
(int.from_bytes(SHA256(f"{domain}/beta/{counter}"),"big") mod
(2^19-1))`. Examine the stream in increasing counter order. Exclude beta 3,
duplicates and any candidate in the Frobenius orbit of beta 3 or an already
admitted candidate. A candidate first passes these **pre-outcome** checks:

1. Rank of all 19 conjugates `beta^(2^j)` is 19, trace is 1, and
   `V_i=span_F2{beta^(2^i),beta^(2^(6+i))}` has rank 2 for each `0<=i<6`;
   their total rank is 12.
2. The four x values in `V_0` lift, under exact group law, to exactly seven
   points `F_0`, as in beta 3. The other five factors are Frobenius rotations
   and must also have seven points.
3. The *single-factor* cofactor-projected column `[4]F_0` has the same number
   of distinct points as beta 3's `[4]F_0` in #767. This is computed only
   from the seven factor points, never a tuple sum or target-support oracle.
   Record both this **signed** point count and the quotient count under
   `P ~ -P` (including O); the signed count, not the quotient, is the
   primary matching condition. This avoids hiding a changed physical column
   behind a sign-paired count.

Primary selection takes the first four candidates passing all three checks.
Stop scanning after the fourth. If fewer than three primary candidates occur
in all 4096 counters, use the saved counter-order preflight stream to append
the first candidates satisfying checks 1–2, without the equal-column check,
that are outside every already selected Frobenius orbit, up to four total.
If fewer than three still exist, preserve the receipt and run only those
admitted. In all cases, archive *every examined counter*, rejection reason,
F0 size, projected-column count when available, and selected list. Source,
selection receipt, #767 reference hash and caps must be committed in this PR
before any candidate six-sum support measurement. The preflight may calculate
normality, factor lifts and `[4]F_0` only. It must not call a six-sum oracle,
read the #767 target verdicts for candidate decisions, or see solver outcomes.

## Exact support and target comparison

For each frozen beta, build six rotated 7-point factor lists exactly as in
#767. Enumerate all `7^6=117649` labelled tuples by group law and save
complete full-point and `[4]`-projected multiplicity histograms, factor
lists and witnesses. Visit *all* q subgroup targets by accumulating the
generator `[4]H` from the identity; save a q-entry little-endian unsigned
32-bit projected multiplicity array in scalar order, including zero entries,
and verify its sum is 117649. No sampling certifies a miss.

Read #767's eight n19 point labels byte-for-byte from its pinned archive.
For each, report projected `R=[4]Q` multiplicity, all four full-point
`Q+T` torsion-coset counts and exact witnesses where present. A planted label
for beta 3 may be negative for another beta; the historical class label
remains unchanged. Do not replace any target. Compare each candidate to the
reference's all-q projected support: exact support and miss counts, labelled
and full/projected duplicate counts, exact projected energy
`E = sum_R multiplicity(R)^2`, unordered colliding tuple-pair count
`(E-7^6)/2`, and Cauchy–Schwarz effective-support floor `(7^6)^2/E`
(as an exact fraction), factor and signed projected-column sizes, intersection,
union, candidate-only and beta3-only targets, common misses, and the fixed
target verdicts. Ratios of support counts are support ratios, never speedups.

An independent bit-serial/Fermat field and separate curve law must rebuild
each factor and **directly enumerate every labelled tuple** (no producer
histogram accumulation), compare both complete multiplicity maps and the
all-q array, check witnesses, and recompute every fixed Q+T count. It must
also independently replay selection from the frozen hash schedule and compare
the preflight receipt, without using tuple-support outcomes for selection.
The parent archive hash and source/input hashes are fail-closed inputs.

## Caps, accounting and archive

Preflight has a 60-second wall and 128-MiB peak-RSS acceptance cap; each
candidate producer has 300 seconds and 512 MiB, and each independent verifier
has 600 seconds and 512 MiB. Process exit, `SIGALRM` wall interruption and
post-run RSS acceptance gates are recorded. Charge field/basis setup,
factor lifting, all tuple and subgroup-target group operations, projection,
reference comparison, target replay, serialization and independent verification.
Report operation types separately as the primary accounting and wall/CPU/RSS
as host-specific diagnostics. Preserve partial files, stdout, stderr and a
failure receipt if any child fails or exceeds a cap; do not substitute an arm.

The source/input/preflight freeze is a dedicated commit and draft PR before
measurement. A second commit records raw archive, SHA ledger, receipt, concise
result, decision and the canonical
`docs/index-calculus-scoreboard.html` stage-only rows. CI performs hash-only
preflight before evidence exists and archive-only independent replay once it
does; it must never run fresh support outcomes. Review the final diff and
exact-head CI, then merge under `AGENTS.md`'s standing workflow.

## Exact pre-outcome source and input freeze

The 25 examined hash counters admitted beta `338435` (counter 1), `303097`
(counter 12), `464276` (counter 13) and `42605` (counter 24). All four pass
normal rank 19, seven physical F0 points, seven signed projected F0 points
and four projected sign-orbits including O. The fallback was not used. The
independent bit-serial preflight replay agrees on every counter, rejection
reason and admitted beta. This selection precedes every candidate six-sum
support outcome. `selection.json` SHA-256 is
`b88d2ca8c33b90650ac466937a01e82d428aa6bd1ef53f6da7f5958c843a5d6d`.

The sweep `FROZEN.json` SHA-256 `44ac5784e3544c3bc32ed5b63f0c93ff27029c6a12b4598a26160cf492516543`.
The runner and CI check that independent literal anchor and then every frozen
source, input, selection, parent-source and #767 archive hash. From the
repository root, the **pre-outcome** check is:

```sh
python3 research/notes/ecc2k130/rotated_beta_sweep_20260925/ci_replay.py
```

After exact-head pre-outcome CI is green, the single frozen support run is:

```sh
python3 research/notes/ecc2k130/rotated_beta_sweep_20260925/run.py --out /private/tmp/rotated-beta-sweep-run-20260925
```

The first command reads only source, selection and pinned reference bytes; it
never computes a new candidate six-sum support outcome. The later archive CI
invocation uses `--evidence research/notes/ecc2k130/rotated_beta_sweep_20260925/evidence`
and independently re-enumerates only the committed outcomes.

# ECDLP index-calculus frontier campaign — protocol

**Status: DRAFT, not frozen.** Freezes when the Phase 0 pilot is committed
and the per-size parameters below are filled in from it. Nothing measured
under Phase 1 or 2 may run before the frozen protocol is committed.

**Purpose.** Replace the patchwork frontier of
`research/ecbench_frontier_20261013/` with one campaign whose numbers are
comparable across methods *and sizes*, for §7 of the survey paper
(`papers/ecdlp-index-calculus-program/`). The question:

> For each method family, what does solving **one previously unseen public
> target** cost, cold and online, as a function of `log₂ r`, and where (if
> anywhere) does an index-calculus arm cross the best rho?

**Class.** Measurement / accounting. No method is changed by this campaign.
Public synthetic curves and hash-to-subgroup targets only; no deployed-curve
or key-recovery claim.

## Why a new campaign

The current frontier (2026-10-13 synthesis) cannot answer the question:

| defect | effect | fix here |
|:--|:--|:--|
| each study used its own curves | ratios not rankable across studies | one fixed curve ladder per family, every arm on every curve |
| one or two sizes per comparison | no exponent established | ≥ 5 sizes per family; exponent fits with intervals |
| arms with unpriced native work | `S` is a lower bound, cannot be ranked against complete counts | Phase 0 pricing audit; unpriced arms reported as bounds with a sensitivity band, never ranked as complete |
| mostly L0 macOS | wall times descriptive | counts on any host; wall only from L2 runs on RunPod |
| cross-host counter drift (inherited F4) unresolved | F4 stage counts not replay-stable | F4 arm runs only after the drift is explained, or is reported per-host |
| tuned on the measured targets in places | selection bias | parameters chosen on a disjoint *train* panel, frozen, then measured on *hold-out* targets |

## Families, curves and arms

Sizes are chosen so that every arm of the family finishes at the smallest
sizes and rho finishes at all of them; an IC arm that exhausts or times out
at a size is a row with that status, never dropped.

### K — Koblitz `E_a : y² + xy = x³ + ax² + 1` over `F_{2^n}`

`n ∈ {13, 17, 19, 23, 29, 31, 37, 41}` (both `a` where the subgroup is
large; the curve id is the ecbench identity, never the bit count).

| arm | method | role |
|:--|:--|:--|
| rho-strong | `rho.signed_frobenius_strong` | reference (incumbent bound `ECBND1h05a2ea918f3c`) |
| rho-frob | `rho.signed_frobenius` | baseline |
| bsgs-neg | `bsgs.negation` | generic baseline |
| claw | `claw.pair_table` | generic baseline with table |
| ic-compact | `ic.pipeline`, `compact-orbit-scan` + `mitm-frobenius-counted:m=2` | IC candidate (best lower-bound arm to date) |
| ic-negfold | `ic.pipeline`, `koblitz-standard-subspace` + `pdp3-koblitz:m=2,engine=rootfind-negfold` | IC candidate (best complete arm to date) |
| ic-f4 | `ic.pipeline`, `pdp3-koblitz:m=3,engine=inherited-f4` | algebraic baseline (gated on the drift fix) |
| ic-shared | `ic.shared_rank` | online / amortised candidate |
| rho-strong-aa | as reference | A/A control |

### P — prime-order curves `y² = x³ − 3x + b` over `F_p`

`bits ∈ {20, 24, 28, 32, 36, 40}`: the deterministic P-256-class rungs of
`prime_fast::find_a3_curve` (largest `p ≡ 3 mod 4` below `2^bits`, smallest
`b` giving prime order, order by Hasse-interval BSGS), listed in
[`prime_curves.jsonl`](prime_curves.jsonl) by
`examples/ic_frontier_prime_curves.rs`. (`prime_search` counts points in
`O(p)` and does not finish at 28 bits; Amendment A1.) The PKM arm needs
2-adicity ≥ 12 and runs on its own registered ladder,
`research/pkm_tower_ecbench_20261013/curves.jsonl` (20–28 bits).

| arm | method | role |
|:--|:--|:--|
| rho-neg | `rho.negation` | reference |
| bsgs-neg | `bsgs.negation` | baseline |
| kangaroo | `kangaroo.vow` | baseline |
| ic-pair | `ic.pipeline`, `prime-abscissa` + `mitm:negation_folded=1` | IC candidate |
| ic-rootfind | `ic.pipeline`, `prime-abscissa` + `rootfind-subtract:m=2` | IC candidate |
| ic-pkm | `ic.pipeline`, `prime-kummer:t=T` + `rootfind-mitm` | IC candidate (curves with 2-adicity ≥ 12 only; separate ladder from `examples/prime_tower_curves.rs`) |
| rho-neg-aa | as reference | A/A control |

### Out of scope here, reported from existing studies with their own boundary

N83 dense root index (not an ecbench method; 8 GiB table), Gaudry `k = 4`
over `F_{p⁴}` (modelled crossover), GPU rho. The paper cites them separately
and does not place them on the campaign's axes.

## Phases

### Phase 0 — pilot and pricing audit (any host, counts only)

1. Build the release `ecbench`; record its SHA-256 and `ecbench methods`.
2. **Pricing audit**: for every IC arm, run one execution per family at the
   smallest size and list every `*_uncharged` counter. For each, record
   whether a pinned conversion exists. Arms with none stay lower bounds;
   for each such counter state the conversion that would complete it and
   a sensitivity band (cost if the counter were priced at 1 GAE each).
3. **Parameter selection** on a *train* panel (target seed `20261015100`,
   4 targets × 1 round per candidate setting): factor-base size / columns /
   dimension per size, chosen to minimise mean cold `S`. Candidates per arm
   are listed in `pilot/candidates.json` before the pilot runs.
4. Freeze: fill the per-size parameter table below, commit, then Phase 1.

### Phase 1 — counted frontier (any host; operation counts are host-independent)

- Hold-out targets: `target_kind: public`, target seed `20261015200`,
  **16 targets per curve × 2 rounds**, 1 warm-up, order `alternate`.
- `isolation_required: L0` (counts only), `timeout_seconds` 3,600 per
  execution; a run that hits it is a `timeout` row.
- Linear algebra `incremental-gauss` (stop when the target scalar is
  pinned) for every IC arm, the scalar verified outside the solver
  (Amendment A2).
- `ecbench verify` replay-all on every session; receipts committed.

### Phase 2 — wall time and second host (RunPod)

- The same specs on a RunPod CPU pod with `--cpus auto`; record the
  isolation level each run earns (a container pod is expected to reach L1,
  possibly L2; L3 needs bare metal and is not expected).
- Wall time is admitted only where `ecbench compare` admits it (≥ 5 pairs at
  the required level, A/A interval inside [0.9, 1.1]).
- `ecbench verify` of the Phase 1 sessions on the RunPod host: the replay
  certificate is the cross-host check. Any counter mismatch is a retained
  failed audit and blocks that arm's counts from the bound.

## Analysis (fixed now)

- **Per curve:** mean `S`, `S / floor`, and paired ratio of totals to the
  reference with two-stage cluster-bootstrap 95% intervals (`ecbench
  compare`).
- **Per family and arm:** fit `log₂(S·√r) = α·log₂ r + c` over all sizes the
  arm completed (≥ 4 required to report `α`), bootstrap over targets within
  size; report `α` and `c` with 95% intervals. Rho's expected `α = 0.5`.
- **Crossover:** for an IC arm with `α_IC < α_ρ` (intervals separated),
  report the extrapolated crossover `log₂ r*` with its interval; otherwise
  report "no crossover supported by the data". Extrapolation is labelled as
  such.
- **Phase breakdown:** share of cold GAE per phase (factor base, oracle
  set-up, relations, linear algebra, verification) per arm and size.
- **Online and amortised (secondary):** from the online window and the
  phase split, `T(q) = B + q·O` for `q ∈ {1, 8, 64}` and the break-even
  `q*` against rho, labelled secondary per the ledger rules.
- Lower-bound arms are plotted hatched and never placed above a complete
  arm in a ranking.

## Predictions (stated before data)

- P1: every IC arm's cold `S / S_ρ` is > 1 at every size.
- P2: Koblitz IC arms have fitted `α ≥ 0.5`; their gap to rho does not
  close with `n` on this ladder.
- P3: prime IC arms' `S / S_ρ` grows with size (slope consistent with the
  `p^{1/6}` of the 2026-10-13 PKM study).
- P4: the online window favours `ic-shared` / `ic-compact` at the larger
  Koblitz sizes (constant-factor, secondary).

A result contradicting P1 or P2 on the cold boundary at ≥ 4 sizes would be
the first *advance* of the programme and triggers an independent replay
before any statement.

## Frozen per-size parameters

_To be filled from the Phase 0 pilot._

## Resources

- Phase 0/1: this repository's container (4 vCPU, 15 GiB) or any Linux host;
  counts do not depend on the host.
- Phase 2: a RunPod CPU pod (needs `RUNPOD_API_KEY` in the environment;
  `ecc2k130/runpod_job.py` is the existing GPU-pod launcher and needs an
  S3 bucket for results, so a CPU variant or an SSH pod is required).

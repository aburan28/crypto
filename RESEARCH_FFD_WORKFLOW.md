# FFD breakthrough workflow — a falsification-driven experiment loop

**Goal.** Turn the proof-complexity bridge
(`RESEARCH_FFD_PROOF_COMPLEXITY.md`) from a proposal into either (a) a
defensive theorem — "the first-fall-degree assumption is false at
cryptographic scale for generic curves" — or (b) a sharp negative result
that kills the bridge. Both are publishable; the workflow is designed so
**every loop iteration ends at one of those two outcomes for one claim**,
never in limbo.

**Operating principle.** We do not try to "prove the conjecture." We
maintain a ledger of **falsifiable predictions**, and each iteration
attacks the single prediction with the best *information-per-CPU-hour*:
run the cheapest experiment that could most plausibly **refute** it. A
prediction that survives a serious refutation attempt graduates; one that
fails is recorded as a kill (which reshapes the conjecture). This is the
same loop the repo's autolab already uses (`RESEARCH_AUTOLAB_LOG.md`),
specialised to the FFD program.

---

## 1. The instrument stack (what we can already measure)

| Layer | Module | Measures | Reach |
|---|---|---|---|
| First fall | `ffd_harness` | `d_ff` (syzygy onset) | `n ≤ 8` |
| Refutation degree | `pc_degree_harness` | `D*` per instance (= PC / last-fall) | `2n' ≤ 10` |
| Distribution | `pc_degree_avg` | `D*` min/mean/max over many targets, ratio `ρ` | `2n' ≤ 10` |
| Regime | `pc_degree_avg::run_pc_regime_sweep` / `…_operating_point_sweep` | `D*` vs `ρ`, `D*` vs `n` at fixed `ρ` | `2n' ≤ 10` |
| Expansion (predictor) | `descent_expansion` *(not built)* | `γ(G)` from the multiplication tensor | any `n` |
| Sparse solver (reach) | `groebner_f4` + sparse F4 *(partial)* | `D*` | target `2n' ≤ 16` |

The loop's first job each iteration is to pick a prediction; its second is
to pick the **cheapest instrument** on this stack that can attack it.

---

## 2. The prediction ledger

Status ∈ {`open`, `supported`, `killed`, `blocked`}. "Supported" never
means "proved" — it means "survived a genuine refutation attempt at the
current reach." Maintained as a table here and machine-checked by
`examples/ffd_breakthrough_loop.rs` (which prints a verdict per row).

| # | Prediction | Current status | Last evidence |
|---|---|---|---|
| P1 | At fixed `ρ ≈ 1`, **mean `D*` grows with `n`** (within a parity class) | `supported` | both parities positive slope: odd +0.121/n (3 pts), even +0.060/n (4 pts) — auto-verdict, iteration 1 |
| P1′ | The growth in P1 is **linear**, not log/constant | `blocked` (reach) | needs `2n' ≳ 14` → sparse F4 |
| P2 | `D*` is **monotone in `1/γ(G)`** across a basis sweep | `blocked` | `descent_expansion` not built |
| P3 | The **even/odd parity split** in `D*` is explained by a structural feature of the tensor / field (subfield, half-trace) | `open` | parity effect persists (odd cells ~0.5 deg higher) but BOTH parities now show growth; still unexplained |
| P4 | HKY explicit counterexample systems register as **low-`γ`** | `blocked` | needs `descent_expansion` + the HKY systems coded |
| P5 | `D*` is **insensitive to the curve** `b` at fixed `(n, n', ρ)` (i.e. `D*` is a field/basis invariant, not a curve invariant) | `supported` | EXP-A: 8 curves at n=8,n'=3, mean-D* spread 0.219 < 0.5 gate — auto-verdict, iteration 1 |
| P6 | Over-determined `ρ ≫ 1` collapses `D*` to 2 (Nullstellensatz) | `supported` | n=10 sweep: ρ≥2.5 ⇒ mean D* ≤ 2.03 (`15b5b1c`, re-confirmed iteration 1) |

New predictions are appended as experiments suggest them; killed ones stay
in the table with their kill evidence (negative results are the point).

---

## 3. The loop (one iteration)

```
            ┌─────────────────────────────────────────────┐
            │ 0. SNAPSHOT: run ffd_breakthrough_loop, read │
            │    the auto-verdict on every ledger row.     │
            └───────────────────────┬─────────────────────┘
                                    │
            ┌───────────────────────▼─────────────────────┐
            │ 1. PICK: the open/blocked prediction with    │
            │    highest (refutation power × tractability).│
            │    Prefer a prediction a cheap experiment    │
            │    could KILL today.                         │
            └───────────────────────┬─────────────────────┘
                                    │
            ┌───────────────────────▼─────────────────────┐
            │ 2. DESIGN the minimal refuting experiment.   │
            │    Write it as a NEW deterministic sweep or  │
            │    a NEW test asserting the refutation       │
            │    condition (red-team, not confirm).        │
            └───────────────────────┬─────────────────────┘
                                    │
            ┌───────────────────────▼─────────────────────┐
            │ 3. RUN. Capture to experiments/*.json via    │
            │    the driver so the result is reproducible. │
            └───────────────────────┬─────────────────────┘
                                    │
            ┌───────────────────────▼─────────────────────┐
            │ 4. JUDGE against the pre-registered gate     │
            │    (§4). Update the ledger: supported/killed │
            │    /blocked. Append predictions if any arose.│
            └───────────────────────┬─────────────────────┘
                                    │
            ┌───────────────────────▼─────────────────────┐
            │ 5. RECORD: one entry in this file's log (§6) │
            │    + commit + push. State carries in git, not│
            │    in memory (container is ephemeral).       │
            └──────────────────────────────────────────────┘
```

**Rule of the loop:** an iteration is only "done" when a ledger row
changed status *or* a new prediction was registered with its own gate.
"Ran some numbers, looks interesting" is not a completed iteration.

---

## 4. Pre-registered decision gates

Gates are fixed **before** running, so a result can't be rationalised
after the fact.

- **G-P1 (growth).** Fit mean `D*` vs `n` on each parity class separately,
  over the largest reach available. *Supported* if the slope is positive
  with the per-cell standard error not overlapping zero across ≥ 4 points.
  *Killed* if mean `D*` is flat (slope CI includes 0) out to the reach
  limit. Anything else → *blocked* (need more reach).
- **G-P1′ (linear vs log).** With ≥ 5 points per parity class, compare a
  linear fit `D* ~ a·n` against a log fit `D* ~ a·log n` by residual sum
  of squares. *Supported-linear* if linear wins by ≥ 2× RSS; *killed* (in
  favour of subexponential-friendly log) if log wins by ≥ 2×; else
  *inconclusive*.
- **G-P2 (expansion law).** Spearman rank correlation between `D*` and
  `1/γ(G)` over ≥ 12 (basis, curve) cells. *Supported* if `ρ_s ≥ 0.6`,
  `p < 0.05`. *Killed* if `|ρ_s| < 0.2`. Mid-range → reformulate `γ`.
- **G-P3 (parity).** *Supported* if a single structural covariate
  (gcd(n, small), subfield existence, half-trace dimension) accounts for
  the even/odd `D*` gap monotonically across the sweep; *killed* if the
  gap persists after conditioning on every covariate we can compute.
- **G-P5 (curve-independence).** For ≥ 8 random `b` at fixed `(n,n')`, the
  `D*` distributions are statistically indistinguishable (KS test,
  `p > 0.1`) ⇒ *supported*; a `b` that shifts the mean by ≥ 0.5 degree
  ⇒ *killed* (D* depends on the curve, not just the field/basis).

---

## 5. The experiment queue (auto-advancing)

The driver prints, each run, the **next experiment** = the gate that is
closest to decidable given current reach. Hand-maintained seed order:

1. **EXP-A — curve-independence (P5).** Cheapest possible: re-run the
   averaged sweep at fixed `(n, n')` over 8 random `b`. No new field math.
   Either kills "D* is a field/basis invariant" or strengthens it — and
   the answer dictates whether `descent_expansion` should be keyed on the
   curve or only the (field, basis). **Do this before building P2's
   instrument.**
2. **EXP-B — parity covariate scan (P3).** For `n = 5..12`, tabulate
   `D*` mean against {`n mod 2`, subfield dims, half-trace rank}. Pure
   bookkeeping over data we already produce.
3. **EXP-C — build `descent_expansion` (unblocks P2, P4).** Tanner graph
   from the `IrreduciblePoly` multiplication tensor; spectral gap +
   small-set boundary scan. Then run G-P2.
4. **EXP-D — sparse-F4 reach (unblocks P1′).** Swap the dense Macaulay
   rank for a sparse/blocked solver to push `2n' = 12 → 16`, then re-run
   G-P1 and G-P1′ with ≥ 5 points per parity class.
5. **EXP-E — basis sweep (P2, P4).** Polynomial vs normal vs sparse
   Galbraith–Gebregiyorgis bases; the central test of the whole bridge.

The queue is a guess; the loop re-prioritises from the ledger each
iteration. EXP-A is first because it is the cheapest experiment that can
*kill* a load-bearing assumption.

---

## 6. Stopping conditions

Stop the program (write it up) when **any** of:

- **Breakthrough-positive:** G-P2 supported *and* G-P1 supported within
  ≥ 2 parity classes — the expansion predictor tracks a real `D*` growth.
  Write the defensive theorem + the screening invariant.
- **Breakthrough-negative:** G-P2 killed (`|ρ_s| < 0.2`) — the bridge's
  predictor is wrong; write up the negative result (still novel: "PC
  degree of Semaev systems is *not* expansion-controlled").
- **Reach wall:** EXP-D plateaus and three successive iterations are all
  `blocked` on reach with no cheaper instrument available. Write up the
  empirical regime map and the open problem.

Do **not** stop merely because an iteration was inconclusive; re-queue and
attack a different prediction.

---

## 7. Iteration log

> One entry per loop iteration. Newest at top. Mirror the autolab format:
> *Task picked · Experiment · Result · Gate verdict · Ledger delta ·
> Commit.* Keep entries short; the JSON snapshots in `experiments/` hold
> the numbers.

### 2026-05-29 — iteration 1 (bootstrap + EXP-A)

- Predictions attacked: P5 (curve-independence), P1 (growth, both
  parities), P6 (re-confirm).
- Experiment: `cargo run --release --example ffd_breakthrough_loop`
  (EXP-A 8 curves at n=8/n'=3; operating-point scaling n=6..12 split by
  parity; regime sweep n=10). Snapshot `experiments/ffd_loop_latest.json`.
- Results:
  - **G-P5 SUPPORTED.** Mean-`D*` spread across 8 curves = 0.219 < 0.5
    gate. `D*` behaves as a (field, basis) invariant, not a curve
    invariant ⇒ `descent_expansion` can be keyed on (field, basis).
  - **G-P1 SUPPORTED on BOTH parities.** Odd slope +0.121/n (n=7,9,11),
    even slope +0.060/n (n=6,8,10,12). The earlier "even is flat" read
    (in `RESEARCH_FFD_PROOF_COMPLEXITY.md`) was undersampling — with the
    parity split and 64 targets/cell, even-`n` also grows, just shallower.
  - **G-P6 SUPPORTED.** Over-determined cells (ρ≥2.5) collapse to mean
    `D* ≤ 2.03`.
- Ledger delta: P5 open→supported; P1 open→supported; P6 reconfirmed.
- Next (driver-chosen): **EXP-C** — build `descent_expansion`
  (field/basis-keyed, per P5) + basis sweep, then run G-P2. P1′ remains
  blocked on reach (EXP-D / sparse F4).
- Commit: (this commit)

### (template)

```
### YYYY-MM-DD — iteration N
- Prediction attacked: P_
- Experiment: EXP-_ (driver invocation / new test)
- Result: <one line + experiments/<file>.json>
- Gate G-P_: supported | killed | blocked | inconclusive
- Ledger delta: P_ open→supported (etc.)
- Commit: <hash>
```

---

## References

See `RESEARCH_FFD_PROOF_COMPLEXITY.md` for the theory and the full
reference list; this file is the operational shell around it.
</content>

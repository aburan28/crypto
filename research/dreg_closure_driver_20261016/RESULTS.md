# Results: the closure's rise is not driven by N; ℓ = 7 or m·ℓ ≥ 21 remain

Scored natively by `dreg_score scale runs/*.d*.jsonl` (`score-output.txt`)
against `PREREGISTRATION.md` (commit `b85dd3226`).

- **Records.** Every value comes from a committed record in `runs/`.
- **Run.** 2026-10-16, 20:49–23:06 UTC. 20 processes ran (16 draws and 4
  controls), all exited 0, and none was killed.
- **Binaries.** These are the scale study's binaries, verified against
  `runs/binary.sha256`.
- **Run report.** With graph and PDF: `reports/run-20261016T2049Z/`.

## The measurement

There are 16 unsatisfiable draws, four per cell, each checked by exact
solution count. Degree 4 was built in full on every draw, and every closure
saturated.

| m | (n, ℓ) | N | m·ℓ | closure degree per draw | cell | XL refutes at 4 | log₂ ops ÷ log₂ cols | peak RSS |
|--:|---|--:|--:|---|---|--:|---|--:|
| 3 | (22, 6) | 40 | 18 | 4 4 4 4 | flat | 0 | 2.28–2.29 | 5.1 GB |
| 3 | (24, 6) | 42 | 18 | 4 4 4 4 | flat | 0 | 2.29 | 7.8 GB |
| 4 | (11, 5) | 42 | 20 | 4 4 4 4 | flat | 0 | 2.18–2.25 | 4.0 GB |
| 4 | (12, 5) | 44 | 20 | 4 4 4 4 | flat | 0 | 2.22–2.33 | 7.8 GB |

- **Q7 holds.** Both m = 3, ℓ = 6 cells refute at 4 at N = 40 and 42, the
  unknown counts where ℓ = 7 did not refute (scale study).
- **Q8 is refuted.** Both m = 4 cells at m·ℓ = 20 refute at 4.
- **Q5 and Q6, carried over, hold.** XL refutes at 4 in 0 of 16 draws, and
  the work ratio is 2.18–2.33.
- **Soundness control.** It replayed four satisfiable draws: m = 4 (11, 5)
  draws 0 and 1, and (12, 5) draws 0 and 1. Neither XL nor the closure
  refuted any of them. The m = 3 cells recorded no satisfiable draw at
  their positive surplus, as registered.

**Verdict, by the registered table: row 1, "ℓ = 6 cells flat; m = 4,
ℓ = 5 cells flat".**

- **N:** rejected. At N = 40 and 42, ℓ = 6 refutes at 4 and ℓ = 7 does
  not. The unknown count does not set the degree.
- **m·ℓ with a threshold ≤ 20:** rejected. m·ℓ = 20 refutes at 4.
- **Surviving:**
  - the degree rises at **ℓ = 7**;
  - the degree rises once **m·ℓ ≥ 21**.
  - Every cell of both studies (56 draws) fits both.
  - At m = 3 they coincide: ℓ = 7 is m·ℓ = 21. They differ at m = 4,
    where m·ℓ ≥ 21 is first reached at ℓ = 6 (m·ℓ = 24).

## What this changes for ECC2K-130 (derived)

- **Both surviving explanations place the m = 4 rise well below the
  priced ℓ.**
  - The "ℓ = 7" explanation puts the first ℓ at which m = 4 exceeds 4 at
    7. The "m·ℓ ≥ 21" explanation puts it at 6.
  - Either way it is far below the ℓ = 33 at which
    `RESEARCH_ECC2K130_DESCENT_DEGREE.md` §6–§7 prices the m = 4
    "4, constant" row (`2^68.43`).
- **So that row loses its support.** It is downgraded from "open
  extrapolation" to **contradicted by both surviving explanations**. Not
  measured at m = 4, ℓ ≥ 6.
  - The degree-5 rows become the optimistic case: m = 4 at `2^71.13`
    (`2^+10.32` over rho), and m = 3 at `2^85.18`.
  - Degree 5 itself is unmeasured above N = 25. No row beats rho, and the
    free-oracle floor still binds.
- **Against the success criterion.** This narrows which extrapolation the
  model may use. It is not an improvement over the baseline.

## Successor

- **The cell that separates the two survivors** is m = 4, ℓ = 6
  (m·ℓ = 24). "m·ℓ" predicts it rises; "ℓ = 7" predicts it holds.
  - Under 155k columns it needs N ≤ 44, so n ≤ 10, surplus ≤ −14.
  - At that surplus, unsatisfiable draws are very rare. The probe here
    found none in 243 draws even at S = −7 for m = 3, ℓ = 7.
  - It needs either a larger memory host or a sparse closure.
- **For the cost rows,** the open question is now degree 5: does it hold
  past ℓ = 7? That also needs a sparse closure (degree 5 at N = 40 is
  760k columns).

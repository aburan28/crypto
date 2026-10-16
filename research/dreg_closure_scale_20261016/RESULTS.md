# Results: the MutantXL closure stays at degree 4 to N = 42 at m = 4, but rises past 4 at m = 3, ℓ = 7

Scored natively by `dreg_score scale runs/*.d*.jsonl` (`score-output.txt`)
against `PREREGISTRATION.md` (commit `f10c86d38`).

- Every value comes from a committed record in `runs/`.
- Run: 2026-10-16 17:05–18:09 UTC, 41 processes (40 draws plus 1
  control), all exit 0, none killed.
- Binaries: `runs/binary.sha256`.
- Run report with graph and PDF: `reports/run-20261016T1705Z/`.

## The measurement

There are 40 unsatisfiable draws, four per cell, each checked by exact
solution count. Degree 4 was built in full on every draw, and every closure
saturated.

| m | (n, ℓ) | N | closure degree per draw | cell | XL refutes at 4 | log₂ ops ÷ log₂ cols | peak RSS |
|--:|---|--:|---|---|--:|---|--:|
| 3 | (11, 5) | 26 | 4 4 4 4 | flat | 0 | 2.19–2.20 | 0.1 GB |
| 3 | (13, 5) | 28 | 4 4 4 4 | flat | 0 | 2.21–2.22 | 0.3 GB |
| 3 | (15, 5) | 30 | 4 4 4 4 | flat | 0 | 2.23–2.24 | 0.5 GB |
| 3 | (16, 6) | 34 | 4 4 4 4 | flat | 0 | 2.26–2.27 | 1.4 GB |
| 3 | (18, 6) | 36 | 4 4 4 4 | flat | 0 | 2.27–2.28 | 2.1 GB |
| 3 | **(19, 7)** | 40 | **>4 >4 >4 >4** | **rises** | 0 | — | 3.7 GB |
| 3 | **(21, 7)** | 42 | **>4 >4 >4 >4** | **rises** | 0 | — | 5.8 GB |
| 4 | (12, 3) | 36 | 4 4 4 4 | flat | 0 | 2.18–2.22 | 2.8 GB |
| 4 | (14, 3) | 40 | 4 4 4 4 | flat | 0 | 2.20–2.25 | 6.8 GB |
| 4 | (13, 4) | 42 | 4 4 4 4 | flat | 0 | 2.31–2.33 | 7.1 GB |

">4" means the degree-4 closure saturated with no fall left unmultiplied,
and `1` was not in its row space. At `(19, 7)`, for example, it reached
rank 94,230 of 102,091 columns after 533,998 products.

- **Study verdict, by the registered rule: "rises at N = 40".**
- **Q4 (flat at 4): refuted.** It fails on all 8 draws of m = 3, ℓ = 7,
  and holds on the other 32.
- **Q5 (XL never refutes at 4): holds,** 0 of 40.
- **Q6 (work ratio in 2.0–2.4): holds,** 2.18–2.33 over the 32 refuting
  draws, rising slowly with N.

**What the cells say beyond the registered rule (descriptive, not decided
on).**

- The rise follows ℓ, not N.
  - At N = 40 and N = 42, m = 4 (ℓ = 3, 4) refutes at 4 on every draw,
    while m = 3 (ℓ = 7) does not refute at 4 on any.
  - At m = 3 the change is sharp: 12 of 12 refute at ℓ = 5 (16 of 16 with the earlier study's `(10, 5)`) and 8 of 8 at
    ℓ = 6, then 0 of 8 at ℓ = 7.
- In summand unknowns `m·ℓ`, the refuting cells reach 18 (m = 3, ℓ = 6)
  and 16 (m = 4, ℓ = 4). The rising cells have 21.
- So "the closure degree is set by `m·ℓ`" is consistent with these ten
  cells. It is an unregistered reading, not a finding. The cells that test
  it are m = 4, ℓ = 5 (`mℓ = 20`) and ℓ = 6 (`mℓ = 24`).

**Soundness control.** As registered, it uses the first two satisfiable
draws of each fresh cell with N ≤ 36. The generator recorded only one:
`m4-12-3` draw 0, with 24 solutions.

- XL and the closure built it to degree 4 and neither refuted it.
- The control is thin by construction, and this is stated, not hidden.
- The earlier study's control (0 of 345 satisfiable draws refuted, at
  degree 5) covers the same instrument.

## What this changes for ECC2K-130 (derived)

- **The "4, constant" assumption.** The constant-degree-4 row
  (`RESEARCH_ECC2K130_DESCENT_DEGREE.md` §6) assumed the closure stays at
  4 out to the priced ℓ ≈ 33–45.
  - **At m = 3 that is now measured false from ℓ = 7.** The m = 3
    degree-4 row is withdrawn as an optimistic case.
  - At m = 4 the assumption has held so far, to ℓ = 4 and N = 42.
- **A "5, constant" row** is added to `descent_degree_cost`. Its output is
  `research/notes/ecc2k130/descent_degree_cost-output-20261016-scale.txt`,
  and every earlier row reproduces unchanged.

| m | degree | w | ℓ | N | total | vs rho |
|--:|---|--:|--:|--:|--:|--:|
| 3 | 4, constant (refuted at ℓ = 7) | 1 | 40 | 251 | 2^82.27 | +21.46 |
| 3 | 5, constant | 1 | 41 | 254 | 2^85.18 | +24.37 |
| 4 | 4, constant (holds to ℓ = 4) | 1 | 33 | 394 | 2^68.43 | +7.62 |
| 4 | 5, constant | 1 | 34 | 398 | 2^71.13 | +10.32 |
| 4 | 5, constant | 2 | 34 | 398 | 2^106.52 | +45.71 |

- **Neither constant is measured at the priced ℓ.** At m = 3 the closure
  is above 4 from ℓ = 7. Whether it stops at 5 is unmeasured: degree 5 at
  N = 40 is 760k columns, past the dense engine.
- **The m = 4 rows rest on four ℓ values** against a priced ℓ = 33–34.
- **No row beats rho.** The free-oracle floor still leaves `2^5.71` per
  target, and units are unconverted, as in §2 of the descent note.
- **Against the success criterion,** this measures the closure past the
  previous range. It does not improve the baseline. It tightens one
  derived row (the m = 3 degree-4 row is withdrawn) and extends the m = 4
  row's measured support from N ≤ 32 to N ≤ 42.

## Successor

The question that moves the derived rows is how fast the closure degree
grows with ℓ: in steps of about 7 in ℓ, or with `m·ℓ`.

- **Next, under the same engine:** m = 4, ℓ = 5 (N ≥ 46, degree 4 at
  179k columns, about 15 GB by the quadratic trend of the peaks above,
  so at this container's limit). This tests the `m·ℓ` reading.
- **Beyond that:** degree 5 at N ≥ 40 needs the sparse engine, and the
  repository's sparse refutation engine has no closure step yet.

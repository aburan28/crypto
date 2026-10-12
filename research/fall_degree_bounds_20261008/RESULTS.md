# Results: fall-degree experiments E1, E1b, E2 (`m = 2`)

Scored on 2026-10-08 by `score.py` (output in `score-output.txt`), from
the committed JSONL in `runs/` and `exp/`, against `PREREGISTRATION.md`
and its Amendment 1. E3, E4 and E5 have not run.

## Verdicts

| hypothesis | as run | Amendment 1 | decided by |
|---|---|---|---|
| H1: the trace relation `L` lies in `span(F)` | falsified, 58/112 | **held, 112/112** | the as-run script used constant `Tr(x_3)`; every as-run failure is a constant mismatch |
| C1: `Tr\|_V ≢ 0` gives `d_ff = 2` | **held, 67/67** | same | — |
| C2: `Tr\|_V ≡ 0` with `Tr(b/x_3²) = 1` is refuted at degree 2 | falsified, 8/22 | **held, 17/17** | same cause as H1 |
| C3: `Tr\|_V ≡ 0` with constant 0 gives no trace fall | consistent | random bases: `d_ff = 3` with no degree-1 relation on 4/4 | subfield bases get `d_ff = 2` from Prop. B instead |
| side: `dim(span(F) ∩ R_{≤1}) = 1` when `Tr\|_V ≢ 0`, `n ≥ 12` | **held, 38/38** | same | — |
| H1b: the EXP-J syzygy is `(L+1)·L = 0` | **held, 69/69** | same | — |
| H2a: Semaev `d_F` ≤ planted-linear control `d_F` | **held** at `n = 8 … 18` | same | — |
| H2b: the gap is nondecreasing in `n` | **held**, but flat at 1 for `n = 10 … 18` | same | no growth seen at toy size |

## The E2 table: medians of `d_F`

| `n` | Semaev sat | Semaev unsat | plain sat | plain unsat | planted sat | planted unsat | semi-regular `D_reg` |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 8 | 3 | 3 | 3 | 3 | 3 | — | 3 |
| 10 | 3 | 3 | 4 | 4 | 4 | 4 | 4 |
| 12 | 3 | 3 | 4 | 4 | 4 | 4 | 4 |
| 14 | 3 | 3 | 4 | 4 | 4 | 4 | 4 |
| 16 | 3 | 3 | 5 | — | 4 | 4 | 5 |
| 18 | 4 | 4 | 5 | 5 | 5 | 5 | 5 |

The Semaev values come from `runs/`, random base. The control cells
have 4 draws each, and the Semaev cells 3–8.

## What the results establish

- **The trace relation and its consequences are now proved and
  confirmed.** Prop. A of the ledger
  (`research/notes/index-calculus/RESEARCH_FALL_DEGREE_BOUNDS.md` §5)
  derives `L = Tr(X_1) + Tr(X_2) + Tr(b/x_3²)` from Kosters–Yeo Cor. 4.11.
  E1 confirms it on every draw.
- **The repository's "bounded-defect lemma" is explained.** Its one
  syzygy is the Boolean shadow of the trace relation. That gives a proof
  of existence. It also removes it as evidence for a lower bound: one
  linear equation is all it carries.
- **Subfield bases are easy at degree 2, by a theorem.** Prop. B: the
  quadratic content of the subfield descent lives in an `n'`-dimensional
  space. Degree 2 then yields `n' − 1 + Tr(b/x_3²)` affine relations in
  `X_1 + X_2`, measured on 39 of 40 draws. The outlier has
  `x_3 ∈ F_{2^{n'}}`.
- **The low `d_F` is not generic.** A random quadratic system of the same
  shape, with the trace fall planted, needs one more degree at every
  `n = 10 … 18`. The plain control tracks the semi-regular `D_reg`, and
  Semaev sits one below it. So whatever keeps `d_F` low is structure
  beyond the trace fall, and it is not yet identified.
- **What it does not establish.** Whether that gap grows. At `n ≤ 18` it
  is a constant 1. Kousidis–Wiemers' F4 bound, `d_F ≤ 5` at `n = 48`
  against a semi-regular 8, says it must grow eventually. E3 tests where.

## Next

- **E5 (native closure engine).** It is the gate for E3 and E4. The
  Sage closure costs 3–4 minutes per draw at `n = 22` and will not reach
  `n = 32`.
- **E2 extension.** Run the planted control at `n = 20, 22`, where
  Semaev reads 4 and the semi-regular reference is 5. This is the first
  place the gap could reach 2 within Sage's reach.
- **P3.** Prove that `span(F) ∩ R_{≤1}` is exactly `⟨L⟩` for a random
  base. E1 supports it on 38 of 38 draws at `n ≥ 12`.

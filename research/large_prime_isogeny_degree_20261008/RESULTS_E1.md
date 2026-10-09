# E1 result: eigenvalue orders and the rational-torsion window

Run 2026-10-08 with `eigenvalue_orders.py` at the registered size, every
prime `ℓ ≤ 2²⁰` for P-192, P-224 and P-256.  18 s in all.  Summary in
`results/e1_summary.json`; the per-`ℓ` records (10 MB per curve) are kept
out of git and pinned by `results/e1_raw_sha256.txt`.

Result class: **toy and structural; no ECDLP speedup; `S` unset.**

## Verdicts against the registered predictions

| prediction | P-192 | P-224 | P-256 | reading |
|:--|:--|:--|:--|:--|
| E1.1 model, ratio in [0.5, 2] at R ∈ {4, 8, 12, 24} | **fail** (0.47 at R = 4) | **fail** (0.40 at R = 4) | pass (0.53, 0.85, 0.95, 0.96) | the eigenvalue model over-predicts small orders; the deficit is at R = 4 and shrinks with R |
| E1.2 `r_min = 2` set equals the small primes of the twist order | pass ({23}) | **fail as written** | **fail as written** | the registered statement forgot the ramified primes (3 for P-224; 3, 5 for P-256), which divide the twist order but have no eigenvalue; restricted to Elkies `ℓ` all three agree exactly |
| E1.3 at least 20 Elkies `ℓ ∈ (61, 2²⁰]` with `r_min ≤ 12`, and one with `r_min ≤ 4` above 10³ | **fail** (13; yes, 10453) | **fail** (9; none) | **fail** (12; yes, 10657 and 318281) | the window is real but about half the registered size |
| E1.4 charged crossover, `c_v = 60` | pass | pass | pass | also at `c_v = 20` and `200` |

Counts of `r_min ≤ R` over all Elkies `ℓ ≤ 2²⁰` against the model:

| R | P-192 | model | P-224 | model | P-256 | model |
|--:|--:|--:|--:|--:|--:|--:|
| 4 | — | — | — | — | 4 | 7.6 |
| 12 | — | — | — | — | 19 | 20.1 |
| 48 | — | — | — | — | 59 | 61.2 |

(P-256 shown; the other curves' counts are in the summary file with the
same shape: the ratio rises from about 0.4 at R = 2–4 to about 0.9 by
R = 12.)

## Reading

- The model `Pr[r(λ) ≤ R] = S_R(ℓ)/(ℓ − 1)` treats `λ` as uniform in `F_ℓ^*`.
  The data say small orders are rarer than that by about half at `R ≤ 4`
  and agree within 10% by `R ≥ 12`.  Since `λ μ = p`, the pair is not two
  independent draws; a corrected model should condition on `ord(p mod ℓ)`.
  This is the registered "failing direction" and is reported before any
  further use of the model.
- E1.2's failure is a wording error in the protocol, not in the data: the
  registered set should have been "Elkies primes dividing the twist order".
  The Elkies-restricted identity holds exactly on all three curves, and
  ties the instrument to the walker's `twist_factors`.
- E1.3 fails on the count but the window exists: P-256 has Elkies degrees
  179, 181, 271, 373, 617, 1289, 2647, 10657, 67891, 76213, 169321, 318281
  with a kernel over `F_{p^r}`, `r ≤ 12`, and at `ℓ = 10657` and `318281`
  the extension degree is at most 4.  Those are the degrees where C4 is
  the cheap route and `Φ_ℓ` is never built.

## Decision, per the registered rule

E1.3 failed, so C4 stays a curiosity for the walk as a whole and the
large-`ℓ` question is decided between C1 and C5.  C4 is still admitted
for the specific degrees listed, since E1.4 holds there, and a walker
edge kind with an `F_{p^r}` kernel-point witness is worth specifying for
them.

# ECDLP index-calculus frontier campaign (2026-10-15)

**Decision: measurement / accounting.** One campaign, every method family on
the same fixed curve ladders, every cost priced across the AGENTS.md §5a band.
**No index-calculus arm beats generic rho on any curve at any band weight; no
index-calculus arm has a fitted cost exponent below rho's.** Operation counts
only (L0); wall time and a second-host replay certificate are Phase 2 (RunPod,
not yet run). Public synthetic curves and hash-to-subgroup targets; no
deployed-curve or key-recovery claim.

| | |
|:--|:--|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) (frozen before Phase 1; amendments A1–A6 dated) |
| pilot (Phase 0) | [`pilot/FINDINGS.md`](pilot/FINDINGS.md), [run report](reports/run-20261015T0436Z/report.pdf) |
| parameter selection | [`select/SELECTION.md`](select/SELECTION.md), [`params.json`](params.json) |
| Phase 1 | specs `phase1/spec-*.json`, sessions `phase1/sessions/`, audits `phase1/audits/` (19/19 ok, every replay identical), [run report](reports/run-20261015T0709Z/report.pdf) |
| analysis | [`analyse.sh`](analyse.sh) (native `ecbench band`, `ecbench bound fit --weight --param-rule`), [`analysis/`](analysis/) |
| paper | `papers/ecdlp-index-calculus-program/` §10.1 |

## Results

Mean `S` of each index-calculus arm over the cheapest constant-memory generic
arm on the same curve (`rho-frob` Koblitz, `rho-neg` prime), at w = 0 / 0.1 / 1
GAE per native operation; 16 public targets × 2 rounds, every run verified
([`analysis/ratios_to_rho.tsv`](analysis/ratios_to_rho.tsv)):

| curve | log₂ r | arm | w=0 | w=0.1 | w=1 |
|:--|--:|:--|--:|--:|--:|
| Koblitz n=13 | 11.0 | compact orbit | 6.9 | 9.3 | 27.1 |
| Koblitz n=17 | 16.0 | compact orbit | 3.8 | 4.9 | 11.8 |
| Koblitz n=19 | 18.0 | compact orbit | 4.9 | 7.5 | 23.6 |
| Koblitz n=23 | 22.0 | compact orbit | 5.5 | 7.4 | 18.4 |
| Koblitz n=37 | 27.8 | compact orbit | 5.2 | 8.3 | 23.3 |
| Koblitz n=43 | 32.1 | compact orbit | 7.4 | 10.2 | 23.1 |
| Koblitz n=47 | 36.6 | compact orbit | 9.9 | 12.7 | 25.3 |
| Koblitz n=41 | 39.0 | compact orbit | 15.0 | 17.8 | 30.3 |
| prime 20–40 bit | 20–40 | pair table | 9.8 → 95.2 | 10.5 → 101.8 | 13.5 → 131.4 |
| PKM 20–28 bit | 18.9–26.2 | Kummer tower | 30.3 → 141.3 | 30.4 → 158.7 | 31.0 → 238.7 |

Negfold root finding: 12–45× rho to n = 23, about 1,900× at n = 37, does not
complete within 5·10⁷ trials from n = 43 (A4b). Inherited F4 (n ≤ 19):
`S` 16 → about 2·10⁴ at n = 13 once its word XORs are priced at w = 0.1.

Fitted exponents `α` (total cost ∝ r^α, 95% intervals, ≥ 4 sizes per tier):
prime pair table 0.620 [0.557, 0.663] vs rho 0.435 [0.331, 0.501] at w = 0
(separated at w = 0 and 0.1, overlapping at w = 1); Koblitz compact orbit
n = 37–47 0.589 [0.518, 0.743] vs signed-Frobenius rho 0.436 [0.358, 0.508]
at w = 0 (overlapping at w = 0.1 and 1); every other IC/rho pair overlaps.

## Limits

L0 counts only; parameters selected on four training targets per size;
`ic-shared` leaves three categories uncounted (lower bound only); the medium
prime tier has two sizes; the band weighs a 64-bit word XOR like a table
lookup, which is pessimistic for F4.

**Correction (2026-10-15).** An earlier version of this README said `ecbench
bound check` of the committed `docs/bounds/records/` failed independently of
this work. It did not: adding `weight` and `param_rule` to `FitOptions`
serialised them as `null` into every re-derived bound, which changed every
bound id. They are now omitted when unset; all 46 committed bounds re-derive
with their recorded ids, and the campaign's 63 fits were refitted (values
unchanged, ids changed).

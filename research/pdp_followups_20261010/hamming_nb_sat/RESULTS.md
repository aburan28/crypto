# Results: Hamming-weight normal-basis factor base under off-the-shelf SAT

Protocol: [`PROTOCOL.md`](PROTOCOL.md), frozen before the run.  Frozen
output: [`results/n*.jsonl`](results/) (one header line and sixteen cell
lines per file; 8 files, **128 cells**), the generated table
[`results/summary.md`](results/summary.md) (`sh summarize.sh`, `jq` over the
JSONL) and `results/run.log`.  Instrument:
[`examples/hamming_nb_sat_panel.rs`](../../../examples/hamming_nb_sat_panel.rs)
as committed in `5e208e7`; the later rustfmt and clippy commits change no
behaviour, checked by replaying both `n = 9` files with the current binary
(identical verdicts, conflict counts, targets and relations).  Solvers:
CryptoMiniSat 5.11.15 (native XOR rows) and CaDiCaL 1.7.3 (Tseitin-split
parity), Ubuntu 24.04 packages.  Host: shared 4-vCPU container, two panel
streams running concurrently; wall time is a practicality note only.

**Stage diagnostic.**  `S`, end-to-end cost and speedup are unset.  Nothing
here is a statement about `n = 131`.

## 1. The table

Per file and solver, over 8 targets (6 uniform, 2 planted): `T_enum` is the
`m − 1`-summand enumeration floor (`|F_w|` for `m = 2`, `C(|F_w|, 2)` for
`m = 3`); "oracle" is the exact pair-table decision (build time excluded).
CryptoMiniSat reports no counters when `timeout` kills it, so its conflict
columns cover decided cells only.  Medians are lower medians over decided
cells.

| n | w | m | `|F_w|` | `T_enum` | oracle, max per target | solver | found | unsat | timeout | oracle: decomposable | conflicts on found | conflicts on unsat | median wall found |
|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|:--|:--|--:|
| 9 | 3 | 2 | 127 | 127 | < 1 ms | CryptoMiniSat | 8 | 0 | 0 | 8 | 50 … 300 | — | 8 ms |
| 9 | 3 | 2 | 127 | 127 | < 1 ms | CaDiCaL | 8 | 0 | 0 | 8 | 48 … 791 | — | 11 ms |
| 9 | 3 | 3 | 127 | 8,001 | < 1 ms | CryptoMiniSat | 8 | 0 | 0 | 8 | 258 … 5,474 | — | 32 ms |
| 9 | 3 | 3 | 127 | 8,001 | < 1 ms | CaDiCaL | 8 | 0 | 0 | 8 | 302 … 8,713 | — | 41 ms |
| 11 | 3 | 2 | 199 | 199 | < 1 ms | CryptoMiniSat | 8 | 0 | 0 | 8 | 237 … 1,589 | — | 18 ms |
| 11 | 3 | 2 | 199 | 199 | < 1 ms | CaDiCaL | 8 | 0 | 0 | 8 | 724 … 12,483 | — | 98 ms |
| 11 | 3 | 3 | 199 | 19,701 | < 1 ms | CryptoMiniSat | 8 | 0 | 0 | 8 | 619 … 280,932 | — | 965 ms |
| 11 | 3 | 3 | 199 | 19,701 | < 1 ms | CaDiCaL | 8 | 0 | 0 | 8 | 434 … 49,458 | — | 604 ms |
| 23 | 2 | 2 | 323 | 323 | 12 ms | CryptoMiniSat | 2 | **6** | 0 | 2 | 14,596 … 19,091 | 43,516 … 56,271 | 1.3 s |
| 23 | 2 | 2 | 323 | 323 | 12 ms | CaDiCaL | 2 | **6** | 0 | 2 | 18,223 … 28,100 | 37,637 … 48,472 | 1.7 s |
| 23 | 2 | 3 | 323 | 52,003 | 1.65 s | CryptoMiniSat | 0 | 0 | **8** | 6 | — | — | — |
| 23 | 2 | 3 | 323 | 52,003 | 1.65 s | CaDiCaL | 0 | 0 | **8** | 6 | — | — | — |
| 23 | 3 | 2 | 2,393 | 2,393 | 87 ms | CryptoMiniSat | 3 | 0 | 5 | 4 | 411,835 … 1,330,918 | — | 82.8 s |
| 23 | 3 | 2 | 2,393 | 2,393 | 87 ms | CaDiCaL | 4 | 0 | 4 | 4 | 220,698 … 913,681 | — | 32.6 s |
| 23 | 3 | 3 | 2,393 | 2,862,028 | 362 ms | CryptoMiniSat | 0 | 0 | **8** | 8 | — | — | — |
| 23 | 3 | 3 | 2,393 | 2,862,028 | 362 ms | CaDiCaL | 0 | 0 | **8** | 8 | — | — | — |

Totals: 75 found (every one re-added to its target in the group; ten
algebraic-only models at `n = 11, m = 3` under CryptoMiniSat were blocked
first), 12 unsat, 41 timeouts; **zero disagreements** with the exact oracle.
Every `n = 23, m = 3` cell timed out although the oracle decides the same
targets in at most 1.65 s (`w = 2`) and 0.36 s (`w = 3`); CaDiCaL spent
1.16 … 1.48 M conflicts per 180 s cell at `n = 23, m = 3`.

Instance sizes at `n = 23`: `m = 2`: 618 AND gates and 69 parity rows
(798 variables / 2,140 clauses for CryptoMiniSat; 3,563 / 13,407 for
CaDiCaL); `m = 3`: 2,823 AND gates and 138 parity rows (3,139 / 8,933 and
10,964 / 40,647).  The Python-era Karatsuba circuit at `n = 23` had 8,029
gates.

## 2. Predictions

- **P1 (correctness): PASS.**  128 cells, zero disagreements, 75 verified
  relations, zero spurious relations admitted.
- **P2 (reach: at `n = 23, w = 3, m = 3`, CryptoMiniSat decides at least one
  uniform target within 180 s): FAIL.**  All 6 uniform and both planted
  targets timed out in both solvers.
- **P3 (enumeration: on every refuted uniform cell at `n = 23`, conflicts
  `≥ 0.5 · T_enum`): PASS, by two orders of magnitude.**  The only refuted
  cells are `w = 2, m = 2` (`T_enum = 323`): CryptoMiniSat 43,516 … 56,271
  conflicts (`135 … 174×`), CaDiCaL 37,637 … 48,472 (`117 … 150×`); 3 … 7 s
  per target against the oracle's 12 ms.  At every other `n = 23` setting
  no refutation was proved at all.
- **P4 (scaling: median found wall at `m = 3, w = 3` grows `≥ 2^{1.5}` per
  field bit from `n = 11` to `n = 23`): not scorable.**  There is no found
  cell at `n = 23, m = 3, w = 3`; the censored lower bound is
  `181 s / 965 ms ≈ 2^{7.5}` over 12 bits, i.e. `≥ 2^{0.63}` per bit, below
  the threshold, so the data neither confirm nor refute it.  At `m = 2,
  w = 3` the found-cell medians grow `2^{0.7}` (CaDiCaL) to `2^{1.0}`
  (CryptoMiniSat) per bit, but those medians are over the cells that finished
  (3 and 4 of 8) and are biased low.

## 3. Reading

1. **The native encoding is correct and compact, and reaches no further at
   `m = 3` than the Python circuit did.**  At `n = 11` both solvers decide
   every `m = 3` target in 0.03 … 15 s; at `n = 23` neither decides any in
   180 s, as the Python route decided none in 25 min (the budgets differ, so
   this is "no better shown", not "equal").  At `m = 2` the encoding does
   decide `w = 3` targets at `n = 23` (3 and 4 of the 4 decomposable ones,
   with median 33 s (CaDiCaL) and 83 s (CryptoMiniSat), range 24 … 162 s), which no earlier SAT figure in the repository
   covers.
2. **Wherever a refutation is proved the solver is an enumerator with a large
   constant.**  `w = 2, m = 2` costs `117 … 174×` the `m − 1` floor in
   conflicts and `≈ 320 … 770×` the oracle's wall; at `|F_w| = 2,393` no
   solver proves one in 180 s.  The weight bound enters as a sequential
   counter and the quadratic descended rows do not propagate until nearly a
   whole summand is assigned, which is the Hamming-ideal thread's "tame
   depth within 1–3 of a full summand" seen from the SAT side.
3. **The base is not where the solver wins.**  The one structural edge it was
   expected to have over Gröbner (cardinality constraints natively) did not
   change the shape: found cost still grows by roughly a factor two per
   field bit (`m = 2`), and `m = 3` is out of reach at 23 bits.

**Class: boundary** for off-the-shelf CDCL on this base at `m = 3`, `n = 23`
within 180 s, by the protocol's registered rule (P1 passes, P2 fails).  The
product law is untouched: the base is Frobenius-stable and storage-free,
which the decomposition note's §6 already prices, and the oracle constant it
gives is worse than the pair table's by `10² … 10⁴` at every rung measured.

## 4. Not covered

- **The PB-XOR arm** (LinPB, Yang–Meel CP 2021) and CryptoMiniSat's BNN input:
  unreachable here (`github.com` answers 403; the packaged CryptoMiniSat has
  no BNN parser).  Cardinality is the sequential counter in both arms.
- Budgets longer than 180 s at `n = 23`; `w ≥ 4`; `n = 83` or `n = 131`.
- Isolation: the two streams shared a 4-vCPU container, so wall times carry
  contention and are not speeds.

## 5. Reproduce

```sh
cargo build --release --example hamming_nb_sat_panel
./target/release/examples/hamming_nb_sat_panel --n 23 --w 2 --m 2 \
    --targets 6 --planted 2 --seed 20261011 --timeout-s 180 \
    --out results/n23-w2-m2.jsonl
sh research/pdp_followups_20261010/hamming_nb_sat/summarize.sh   # table from results/
```

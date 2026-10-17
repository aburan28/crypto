# P-256 trapdoor search, round 2 — run report

- **Run id.** run-20261012T2150Z (UTC start 2026-10-12T21:50Z)
- **Source.** `RESEARCH_P256_BACKDOOR_SPECULATION.md` §4, ranks 1–3
- **Code.** `src/cryptanalysis/p256_trapdoor_round2.rs` at this PR's commit;
  `cargo run --release --example p256_trapdoor_round2` (Markdown) and
  `-- --json` (records committed as `records.json`).
- **Host.** Local darwin/arm64 laptop build, `cargo build --release`
  (23 min cold); probe runtime after build: seconds. No isolation level
  claimed: every probe is exact integer arithmetic, deterministic, no
  wall-clock figure is cited.
- **Command.** `cargo run --release --example p256_trapdoor_round2 -- --json`
  with `pollard_iters = 2^22` (reaches prime factors up to ~2^44; `p − 1`
  came out fully factored, so all orders are exact).
- **Runs requested / completed.** 1 / 1, exit 0, no anomaly flags.

## What ran (measured values)

| # | Hypothesis | Probe | Measured |
|---|---|---|---|
| 1 | H1 coefficients | `ord(b) = ord(r) = ord(−3) mod p` | 253 bits, exact; contains the 160-bit prime `Q` of `p−1` |
| 2 | H1 | `ord(j) mod p` | 252 bits, exact; contains `Q` |
| 3 | H1 | residuosity | `b, r` QR non-cubes; `j` QR cube; `−b` flips the QR bit |
| 4 | H1 | Hamming weights | `b` z=+0.81, `r` z=−0.31, `j` z=−0.12, `n` z=+4.88 (n order-derived) |
| 5 | H1 | perfect powers / nearest `c·2^k` | none; gaps ≥ 2^49 |
| 6 | H1 | continued fractions (64 terms) | max partial quotients 1275, 62, 107 |
| 7 | H2 prime | `p−1` factorization | `2·3·5²·11·101·281·1531·65537·490719·6602113·Q₁₆₀` (complete) |
| 8 | H2 | `p+1` factorization | `2^96·7·273921·` 140-bit composite (partial) |
| 9 | H2 | `ord_p(2)` | `(p−1)/2` exactly (2 is a QR; p ≡ 7 mod 8) |
| 10 | H2 | twist embedding degrees (≥2^20 subgroups) | none for k ≤ 100 |
| 11 | H3 endo/isog | GLS gap `4p−t²` | 258 bits |
| 12 | H3 | rational roots of `x³−3x+b` | 1 (⇒ 1 rational 2-isogeny) |

Every value above is copied from `records.json`; nothing is computed in
this report. Anomaly flags: none (threshold z = ±5; the largest
observed is `n` at z = +4.88, recorded but attributed to the `FF`-run
structure `n` shares with `p`, and `n` is point-count-derived, not
seed-derived).

## Graph

`graph.svg` plots the recorded Hamming-weight z-scores of the four
256-bit constants against the binomial expectation (0) and the ±5σ
anomaly threshold, and `factors.svg` the recorded factor-size profiles
of `p−1` and `p+1`. Both plot recorded values only.

## Observations (no interpretation beyond the note)

- No probe flagged an anomaly; see `RESEARCH_P256_BACKDOOR_SPECULATION.md` §5.
- `p−1` was fully factored within the rho budget; `p+1` was not
  (140-bit composite residue) — recorded as partial, not as smooth.

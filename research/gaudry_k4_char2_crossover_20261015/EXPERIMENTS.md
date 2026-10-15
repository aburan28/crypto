# Experiment log: `k = 4` in characteristic 2 past one word (sparse LA, `m = 15–18`, large primes)

Every run made for this study, in order: readings, pilots, failures, claim
runs and audits. It is append-only and updated in the same PR as each run
lands. Times are UTC, 2026-10-15 unless marked. Host: Claude Code cloud
container, Xeon @ 2.10 GHz, 4 vCPU.

Statuses:
- **reading**: existing records, no run;
- **pilot**: before `PROTOCOL.md`; chooses options and predictions, no claim;
- **claim**: registered in `PROTOCOL.md` before it ran;
- **audit**: `ecbench verify --replay-all`.

| # | when | status | what | command / spec | outcome | data |
|--:|:--|:--|:--|:--|:--|:--|
| W1 | 14:50 | reading | Linear-algebra cost in the committed char-2 sessions (#1883, B1–B3), Frobenius arm | `records.jsonl` `linear_algebra` phase | `row_ops` 3.6·10², 8.5–10.7·10³, 4.7–5.3·10⁵, 2.8–3.3·10⁷, 1.9·10⁹ at 14, 62, 262, 1,045, 4,123 columns: ×61 and ×64 per 4× columns, **cubic**. The dense matrix's share of the IC cost: 0.07 % (`m = 12`), 1.06 % (`m = 14`). Projected from `m = 14`, dense LA alone reaches rho's cost near `m = 18`, so the extrapolated crossover of #1883 cannot happen with this matrix | — |
| W2 | ~15:10 | — | `wiedemann` relation matrix added (`ic_framework/linalg.rs`): singleton-filtered core, Wiedemann on `RᵀDR`, Berlekamp–Massey, exact check | `cargo test --release --lib ic_framework::linalg` | 5 passed (planted value pinned on 20 seeds, never before the dense matrix could, at most 16 rows after) | — |
| W3 | ~15:30 | — | One- and two-large-prime collection in `gaudry-d4-bin` (`large_primes`, `small_permille`) | build | compiles; measured in the claim runs | — |
| W4 | ~15:35 | finding | `f2m4_torsion` refuses `m > 15`: group orders, scalars and point keys are one word (`#E ≈ 2^{4m}`). The projected crossover (relation part ÷ 4 per 2 steps of `m` from 6.06 at `m = 14`) is at `m ≈ 17` | `ecbench plan` on `m = 16, 17, 18` | Reported to the user, who chose to widen the harness to 128-bit orders | — |
| W5 | ~16:00–17:30 | — | Widening: `gaudry_quartic_bin` to `m ≤ 20` (`F_q` elements `u32`, orders `u128`, keys mixed to 64 bits past `4m + 1 = 64`, Miller–Rabin for `u128`); `ic_wide.rs`, scalar-generic ports of `rho.negation` and the walk relation loop with `wiedemann`; ecbench routes `f2m4_torsion` with `r ≥ 2^64` there | build; equivalence tests at `S = u64` against the one-word code | see W6 | — |
| W6 | 16:10–16:27 | — | Equivalence tests (`ic_wide`), first run | `cargo test --release --lib -- ic_wide …` | At `S = u64` both ports equal the one-word code; at `S = u128` they **did not** (`walk_operations` 610 vs 498; `total_gae` 8.96·10⁶ vs 6.92·10⁶): `gen_range` over `u128` draws differently. Fixed by drawing through `u64` while the bound fits (`gen_below`) | — |
| W7 | 16:20–16:27 | — | Equivalence tests, after the fix; the touched modules' suites | same, plus `gaudry_quartic_bin`, `ecbench`, `linalg` | **2/2** `ic_wide` (rho at `m = 6, 7, 8`, 3 seeds each; pipeline at `m = 6` with large primes off and with two large primes over half the base, at `u64` and `u128`: identical total, trials, LA work and counters) and 88 others passed | — |
| W8 | 16:28– | audit | Re-audit of #1883's sessions B2 and B3 with this branch's binary (`50f35402…`): do the widened field and orders change any count? | `ecbench verify --replay-all` | running | `audits/reaudit-*` |
| W9 | 16:28– | pilot | Dense against `wiedemann`, and large-prime settings, `m = 10`, 2 targets | `pilot/SPEC-pilot-m10.json`, `--cpus 1 --lock /tmp/bench-a.lock` | dense arm `S` = 71.427 on the first target, **equal to B2's committed value**; `wiedemann` 72.137; one large prime over half the base 217.0 (running) | `pilot/session-m10/` |
| W10 | 16:30–16:45 | failure | Curve search at `m = 16`: `ecbench plan` ran 10 min in `Fq4::new` | `ecbench plan` | The quartic search tries `c₃ = 0` first; those quartics (affine linearised) were all reducible and the scan is `q²` long. Fixed for `m ≥ 16` only (try `c₃ = 1` first), so every `m ≤ 15` field is the one committed. Two `pkill` commands matched their own shell and killed it (no run affected) | — |

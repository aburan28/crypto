# Protocol: the Dickson torus tower, and a deeper-2-adic Kummer control

Lane `dickson-torus`, 2026-10-15. Leads 3 and 4 of
[`docs/ic/factor-base-encyclopedia/README.md`](../../docs/ic/factor-base-encyclopedia/README.md)
§8, with a feasibility derivation for lead 1 (no runs). Background:
[`RESEARCH_PKM_TOWER_ORACLE.md`](../notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md)
(the *PKM note*) §2.1 (towers), §3.4 (degeneracy at infinity), §10–§16
(pilot and rounds 2–7), §10.10 (2-adic table).

This file is committed before any measured cell runs. The runs it lists, and
nothing else, make the results; a change after the first run is a dated,
labelled amendment and both versions are reported.

## 1. Questions

- **Q-A (lead 3).** The Dickson tower `V = {w + w^{−1} : w ∈ u·μ_{2^t}}`,
  step `y ↦ y² − 2`, has two cases (§2.1): *split* (`u ∈ F_p^*`,
  `2^t | p − 1`) and *torus* (`u` in the norm-1 torus `T ⊂ F_{p²}^*`,
  `2^t | p + 1`). Only the split case was ever measured, where it is the same
  F4 computation as Kummer (§10.2). At a fixed `(m, N)`, does the torus case
  give the same F4 solving degree `D` as the split case?
- **Q-A∞.** Is the torus tower degenerate at infinity, like the split Kummer
  and Dickson towers (§3.4)?
- **Q-B (lead 4).** Does deeper 2-adicity `v₂(p − 1)` change the Kummer `D`
  at a fixed `(m, N)`, compared with `p₀` (`v₂ = 18`) and `p₁` (`v₂ = 27`)?
- **Q-C (lead 1, derived only).** What memory and time would `m = 3, N = 21`
  and `m = 4, N = 20` need, projected from the recorded growth ratios?

## 2. Instrument

- **Builder.** `examples/pkm_tower_pilot.rs` with two additions, everything
  else unchanged (presentations, `S₃` chain, targets, seeds, engine, row
  format):
  - `--kinds dickson-torus`: the torus case. `F_{p²} = F_p(i)`, `i² = −1`
    (`p ≡ 3 mod 4`). `ζ = h^{(p+1)/2^t}` for a non-square `h ∈ T`; `u ∈ T`,
    `u ∉ F_p`, `u^{2^{t+1}} ≠ 1` (that is, `u² ∉ μ_{2^t}`);
    `V = {2·re(w)}`, `c = u^{2^t} + u^{−2^t}`, scale 1, steps `D_2`.
    `check_torus` panics unless: `ζ` has exact order `2^t`; `u ∈ T \ F_p`;
    every `w + w^{−1}` (by the general `F_{p²}` inverse) has zero `i`-part and
    equals `V[k]`; every `y² − 4` (`y ∈ V`) is a non-residue (the torus case,
    not the split one); every level `y_j ≠ 0` on `V` (`c` avoids `−2` along
    the chain); `c ≠ ±2`. `Tower::check` then checks `|V| = 2^t` distinct
    values and every tower equation on every element. The kind is appended
    last to the `Kind` enum, so the other kinds' cell seeds are unchanged.
  - `--null-chain`: the shape-matched null at `m ≥ 3`, every `S₃` link of
    the chain given random coefficients on its own support (planted zero
    kept). At `m = 2` the null's draws are unchanged.
  - Unit tests (`cargo test --release --example pkm_tower_pilot`):
    `|V| = 2^t`, `V ⊂ F_p`, `D_{2^t}(y) = c`, `y² − 4` non-residue, `t = 1…12`
    at both torus primes; a split prime is refused; the engine's square
    rules vanish on the torus tower; the chain null keeps its planted zero.
- **Engine.** `f4_fp_tower` (`--engine tower`), the round-2 to round-7
  instrument, reduced presentation, degree bound `n + d + 6`, round-6
  full-rank exit on (default; it changes no `D`, §15). Flags as round 2:
  `--engine tower --budget 7200 --max-nnz 2000000000 --ladder-t none
  --planted 0 --random 2`, seed `0x504B4D54`.
- **Overflow.** The pilot's field arithmetic is `u128`-backed `u64`; the
  engine (and `f4_fp`) assert `p < 2³²`. All primes below are under `2³²`.
- **Checker.** `examples/pkm_tower_check.rs verify` recounts every finished
  tower row by exhaustive search over `V^m`, with `V` rebuilt independently
  from the row's `u, ζ` in `F_{p²}` (new `dickson-torus` branch). `analyze`
  tabulates.
- **Resources.** One thread (`RAYON_NUM_THREADS=1`), `nice -n 10`,
  `ulimit -v 2621440` (2.5 GB address space), one heavy run at a time. A
  system that runs out of memory, time or `--max-nnz` reports no `D` (at most
  a lower bound) and is recorded as out of budget.

## 3. Primes

| name | `p` | form | `v₂(p − 1)` | `v₂(p + 1)` | role |
|:--|--:|:--|--:|--:|:--|
| `T24` | 14,680,063 | `7·2^21 − 1` (23.8 bits) | 1 | 21 | torus, "about 20 bits": the smallest prime `c·2^k − 1` with `k ≥ 20` (`2^20 − 1`, `2^21 − 1`, `3·2^20 − 1`, `5·2^20 − 1`, `7·2^20 − 1` are composite) |
| `T31` | 2,046,820,351 | `61·2^25 − 1` (30.9 bits) | 1 | 25 | torus, 31 bits, under `2³¹` (lazy kernel, as `p₁`); `k = 25` is the deepest for primes `c·2^k − 1` in `(2³⁰, 2³¹)` |
| `p₁` | 2,013,265,921 | `15·2^27 + 1` | 27 | 1 | split Dickson arm (this build), recorded Kummer |
| `p₂` | 3,221,225,473 | `3·2^30 + 1` | 30 | 1 | deep-2-adic Kummer |
| `p₀` | 786,433 | `3·2^18 + 1` | 18 | 1 | recorded pilot prime; positive control |

All five are prime (`factor`). **Lead 4 asked for `v₂(p − 1) ≥ 40`** (for
example `c·2^40 + 1 < 2^62`). Both F4 engines refuse `p ≥ 2³²`
(`f4_fp_tower::Fp::new`, `f4_fp::Barrett::new`: coefficients are stored in
32-bit words), and the largest `v₂(p − 1)` of a prime below `2³²` is 30
(`2^31 + 1 = 3·715827883`). So `p₂` is "whatever the arithmetic safely
supports"; `v₂ ≥ 40` is **blocked by an implementation limit** (engine width),
not attempted, and stays open.

## 4. Cells and runs (in this order)

`m = 2` is `full` (`S₃` plus towers), `m ≥ 3` is `chain`. Two random targets
per system cell unless stated. Projected peak memory is the recorded Kummer
value at the same `(m, N)` (requirements README §3.2; PKM note §11.7): `m = 2,
N = 20`: 473 MB; `m = 3, N = 15`: 1.1 GB; `m = 4, N = 12`: 171 MB. Every
requested cell is under the 2.5 GB cap, so none is skipped for memory.
Projected time, single thread, from the calibration below (about
`2.2·10⁹` multiply-adds per second): `m = 2, N = 20` ≈ 220 s per system;
`m = 3, N = 15` ≈ 1,450 s per system; everything else under a minute.

| run | arm | `p` | `m` | `N` | control | flags beyond §2 |
|:--|:--|:--|--:|:--|:--|:--|
| R0 | positive control: isogeny | `p₀` | 2 | 2…10 | tower | `--kinds isogeny --t-min 1 --max-t 5` |
| R1 | torus | `T31` | 2 | 2…20 | tower | `--kinds dickson-torus --t-min 1 --max-t 10` |
| R2 | torus | `T24` | 2 | 2…20 | tower | same |
| R3 | split Dickson | `p₁` | 2 | 2…20 | tower | `--kinds dickson --t-min 1 --max-t 10` |
| R4 | torus null | `T31` | 2 | 2…20 | null | `--controls null` |
| R5 | torus null | `T24` | 2 | 2…20 | null | `--controls null` |
| R6 | torus | `T31` | 3; 4 | 6, 9, 12; 8, 12 | tower | `--m 3,4 --t-min 2 --max-t-m3 4 --max-t-m4 3` |
| R7 | torus | `T24` | 3; 4 | same | tower | same |
| R8 | split Dickson | `p₁` | 3; 4 | same | tower | same |
| R9 | torus chain null | `T31` | 3; 4 | same | null | `--null-chain --controls null` |
| R10 | torus chain null | `T24` | 3; 4 | same | null | same |
| R11 | deep Kummer | `p₂` | 2; 3 | 2…18; 6, 9, 12 | tower | `--kinds kummer`, `t ≤ 9` at `m = 2`, `t ≤ 4` at `m = 3` |
| R11n | deep Kummer null | `p₂` | 2 | 2…18 | null | |
| R12 | torus | `T31` | 3 | 15 | tower | `--t-min 5 --max-t-m3 5 --targets 0` |
| R13 | torus | `T24` | 3 | 15 | tower | same |
| R14 | split Dickson | `p₁` | 3 | 15 | tower | same, if time allows |
| R15 | torus, second target | `T31`, `T24` | 3 | 15 | tower | `--targets 1`, if time allows |

- **Scaled down, stated now.** At `m = 3, N = 15` one target per prime (R12,
  R13), because a system costs about 25 single-thread minutes; R14–R15 run
  only if the time rule allows. The `m = 3, N = 15` nulls are not run.
  `N = 15` at `m = 3` for the deep-2-adic Kummer is not asked and not run.
- **Time rule.** No run starts later than 150 minutes after R0 started. A run
  not started is reported as not run, never as a result.
- After the runs: `pkm_tower_check verify` on every row file; every finished
  tower row must agree.

## 5. Comparators (recorded, same engine)

Split `D` by `(m, N)`, from committed rows of
`research/pkm_tower_round2_20260925/runs/`: `XV1` (Kummer and split Dickson,
`p₀`, `m = 2`, `N = 2…18`; planted targets up to 14, random above), `XV2`
(`m = 3`, `p₀`), `XV3` (`m = 4`, `p₀`), `K1`/`K1n` (Kummer, `p₁`,
`N = 12…22`), `M3`, `M4` (Kummer `p₁`), `K0` (`p₀`, `N = 20`), `XV5` (Kummer,
`p₂`, planted, `N = 2…16`).

| `m` | `N` | split `D` (Kummer = Dickson where both ran) |
|--:|:--|:--|
| 2 | 2, 4, …, 20 | 3, 4, 4, 4, 4, 5, 5, 5, 5, 6 |
| 3 | 6, 9, 12, 15 | 5, 6, 6, 7 |
| 4 | 8, 12 | 6, 7 |

R3 and R8 re-measure the split Dickson arm on this build at `p₁` (a prime
the split Dickson case never ran at); where they ran, they are the reference.
`m = 3, N = 15` uses R14 if it ran, else the recorded Kummer `D = 7`.

**Calibration disclosed (made before this file, on this build):** Kummer
`p₁`, `m = 3, N = 12` and `m = 2, N = 18`, two random targets each, recorded
cells: `D = 6` and 5, widths 5,019 and 8,644, as recorded; 5.0 s and 19.3 s
per system single-threaded. Smoke: torus `T24`, `m = 2`, `N = 2, 4, 6`, two
random targets: `D = 3, 4, 4`, widths 4, 16, 56. These smoke rows are
scratch and are not results; R2 re-measures them.

## 6. Hypotheses and decision rules (fixed now)

**Cell value.** A run gives a cell `(arm, m, N)` a value when every target it
measured finished (refuted or solved, not timed out or stopped for size) with
the same `D`. Targets disagreeing make the cell *mixed*.

**Q-A.**
- **H_same:** at every `(m, N)` where a torus cell (either prime) and the
  split reference both have a value, they are equal.
- **H_diff:** at least one such pair differs. The differing cells are listed
  with both values.
- Mixed cells are listed and counted toward neither; if any exist the verdict
  reads "H_same (or H_diff) on the valued cells; mixed at …".
- Width and step count equality is reported beside `D`, descriptively ("the
  same F4 computation" only if `D`, width and steps all match).
- **Invalid** if the positive control fails or `verify` contradicts any row.
- The torus `m = 2` line is also read under A1 (§10.8, §11.6): H0, H1a or
  inconclusive. H1a needs a final plateau `L > 10`, which `N ≤ 20` cannot
  give after a rise at `N ≤ 12`; that is stated in advance.
- **Null:** torus `D` against the shape-matched null at each `(m, N)`,
  reported (H1b), gates nothing.

**Positive control (R0).** The instrument separates tower families with
different `D` at the same `(m, N)`: the recorded isogeny `D = 5` at
`m = 2, N = 4…10`, where the split towers give 4. R0 passes if it gives
`D = 5` at those `N` on both targets. If it fails, Q-A and Q-B read invalid.

**Q-B.**
- **H2_same:** Kummer `D` at `p₂` (R11) equals the recorded Kummer `D` at
  `p₀` and `p₁` at every `(m, N)` where both have a value.
- **H2_diff:** at least one differs (listed).
- Scope stated now: `v₂ = 30` against 18 and 27, not P-224's 96; `t ≤ 9`, so
  the 2-adicity beyond `t` is 21 levels at `p₂` (9 at `p₀`).

**Q-A∞ (derived, not measured).** The torus tower's equations are
`y_{j+1} − y_j² + 2` and `y_{t−1}² − 2 − c`: the split Dickson polynomials
with another constant `c`. Their top-degree parts are `−y_j²` and `y_{t−1}²`,
the same as Kummer and split Dickson, so §3.4's argument applies unchanged:
the top-part ideal of the raw system is `(y_{i,j}²)` whenever the `S_{m+1}`
top coefficient `c(x_R)` is non-zero, its Hilbert series is `(1 + t)^N`, and
the Bardet–Faugère–Salvy degree of regularity is `N + 1`. Prediction (derived
from the same fact): since the systems differ from the split Dickson ones only
in constants, H_same is expected; the round tests it.

**Q-C.** Derived only, from the recorded tables, labelled as such: the
per-step growth ratios of width, peak memory, multiply-adds and time at the
last recorded step and as a geometric mean over the recorded range, applied
once more (`m = 3`: 18 → 21; `m = 4`: 16 → 20).

**Not admissible:** changing cells, flags, primes, seeds or rules after the
first run; dropping a target; re-running a system for a different answer.

---
title: "IC next-10, item 3: the PKM tower's solving degree past the memory wall"
date: 2026-10-14
---

# Item 3: the PKM tower's solving degree past the memory wall

**Verdict: INCONCLUSIVE on the plateau, because the wall stopped every run at
`N ≥ 15`. POSITIVE on the conflict with `d + δ`, which is now resolved.**

- **The plateau.** msolve 0.6.5 reproduces the committed F4 degree on every
  control it finished: 15 of 15 systems (`m = 2`, `N = 12–18`; `m = 3`,
  `N = 9, 12`; `m = 4`, `N = 8, 12`), with the same verdict on each.
- At `m = 3`, `N = 21` and `24` it gets only as far as a **productive
  degree-7 step** (`D ≥ 7`). The jobs then hit this host's 45-minute wall
  inside the second degree-7 step, and **no step above degree 7 started**.
  - The design's kill test needs `D ≥ 8` at `N = 21` and `≥ 9` at `N = 24`.
    These runs cannot evaluate it, and they do not show `D = 7` either.
- **The wall is wall-clock time, not memory.** The host ran at a load of
  22–25 on 4 cores, so each single-threaded job got about 19% of a core. The
  peak resident memory was 2.6 GB at `N = 21`, well under the 5 GB cap.
- **The `d + δ` conflict is resolved.**
  - For the `m = 2` full presentation, Proposition 2 is an equality. `δ = N`
    exactly at `N = 12–26`: the inverse `1/S̄` has every one of its `2^N`
    square-free coefficients non-zero. So the Macaulay/XL degree is
    `D_XL = N + 4` (16 → 30), while F4 needs 5–6.
  - A direct Macaulay rank computation confirms `D_XL = 16` at `N = 12`. It
    also gives `D_XL = 13` at `m = 3`, `N = 9`, against F4's 6; `D_XL ≥ 14` at
    `m = 3`, `N = 12`, against 6; and `D_XL = 12` at `m = 4`, `N = 8`,
    against 6.
  - The bound is real and tight, but it bounds XL only. F4 is not bound by it,
    because it reuses rows that reduce to a lower degree.
  - The Hilbert-series predictors fail as the note's §3.4 says they must:
    - the top-part ideal of the full presentation has series `(1+t)^N`, so
      `d_reg = N + 1`;
    - the chain presentation has no finite `d_reg`;
    - the affine Hilbert function of the homogenised generators reaches 0 only
      at `D_XL`.

  None of them predicts the measured `D` at any `N ≤ 18`.
- **A new structural observation.** msolve's F4 trace is identical round for
  round at `N = 18`, 21 and 24 (`m = 3`) through the first degree-7 step. The
  selected pairs, rows, new elements and zero reductions all match. The
  columns match too, except in that step at `N = 18`, and at `N = 21` and 24
  they match there as well. That step adds 2,354 elements at every one
  of these sizes, and `f4_fp_tower` added 2,354 at its step 19 at `N = 18`.
  - So F4's work up to and including its first degree-7 step does not depend
    on `t` once `t ≥ 6`. This is the "level-window" mechanism of the note's
    §10.5, now seen at `m = 3` and in a second engine.
  - What grows with `N` is the degree-7 phase after that step. It has 2 steps
    at `N = 15` and 6 at `N = 18` (`f4_fp_tower`). Whether it ever needs
    degree 8 is exactly what the wall hid.

## 1. Setup

- **Systems.** These are the Kummer tower systems at
  `p₁ = 2013265921 = 15·2^27 + 1`, chain presentation for `m ≥ 3`, in the raw
  presentation (as written, not reduced modulo the tower).
  - The committed systems are from
    `research/magma_f4_comparison_20261013/scripts/`.
  - The new sizes were drawn by `examples/pkm_tower_pilot.rs` with round 2/7's
    flags (`--random 2 --planted 0`, same seed). The run was stopped at once
    with `--max-nnz 1000`, and the systems were exported with
    `--export-magma-rows`.
  - **Identity check:** the draw at `N = 18` (`m = 3`) and `N = 16` (`m = 4`)
    gives the committed targets. The polynomials are byte-identical to the
    committed exports; only the comment line with the row's own verdict
    differs.
- **Conversion and analysis.** `examples/ic_next10_03_pkm_degree.rs` (new;
  Rust) does four jobs:
  - `convert` turns Magma into msolve's format, keeping the variable order;
  - `delta` computes `δ` by a 2-D number-theoretic transform on
    `(Z/2^t)^2`;
  - `macaulay` computes the exact XL degree and affine Hilbert function in the
    tower quotient. No row is ever reused at a lower degree;
  - `top` computes the top-degree parts modulo `y²`, and the Hilbert function
    of that ideal when it is monomial.
- **msolve.** Every run used `msolve -t 1 -v 2`, `nice`, `prlimit --as=5e9`
  and `timeout 2700`.
  - `run.sh` samples peak `VmHWM`/`VmPeak` from `/proc` every second.
  - `parse.sh` reads the `-v 2` round table. `D_msolve` is the largest round
    degree. `D_last_productive` is the largest degree of a completed round
    that added basis elements; this is the lower bound for a stopped run (note
    §16.1).

## 2. Results

### 2.1 msolve against the committed F4 degree

| m | N | targets | committed `D_f4` (`f4_fp_tower`) | msolve | verdict match | peak RSS | wall (s) |
|--:|--:|--:|--:|:--|:--|--:|--:|
| 2 | 12 | 2 | 5 | 5, 5 | yes (`[-1]`) | 6 MB | 1 |
| 2 | 14 | 2 | 5 | 5, 5 | yes | 6 MB | 1 |
| 2 | 16 | 2 | 5 | 5, 5 | yes | 43 MB | 6–7 |
| 2 | 18 | 2 | 5 | 5, 5 | yes | 171 MB | 45–116 |
| 2 | 20 | 1 | 6 | stopped at 45 min in a degree-6 step (`D ≥ 5` productive; degree 6 reached) | — | 681 MB | 2701 |
| 3 | 9 | 2 | 6 | 6, 6 | yes | 6 MB | 1 |
| 3 | 12 | 2 | 6 | 6, 6 | yes | 100 MB | 36 |
| 3 | 15 | 1 | 7 | stopped at 45 min in the 2nd degree-7 step: `D ≥ 7` | — | 1.47 GB | 2701 |
| 4 | 8 | 2 | 6 | 6, 6 | yes | 6 MB | 2–3 |
| 4 | 12 | 2 | 7 | 7, 7 | yes | 477 MB | 595–1242 |

Every finished control agrees in `D` and in verdict: 14 systems, plus the
`m = 3`, `N = 9` pilot system, 15 in all. msolve's sequence of round degrees
has the same shape as `f4_fp_tower`'s: it climbs to the plateau degree, then
falls and finds `1` at a lower degree.

**Instrument check at `p₁`.** The note warns that msolve 0.6.5 failed a
known-answer battery above `2^16` in another project. On these systems, at
`p₁ ≈ 2^31`, every verdict agrees with the committed rows, and those rows were
verified against exhaustive search. The linear-algebra variants `-l 2`
(default), `-l 42` and `-l 44` give identical round tables at `m = 3`,
`N = 12` (`runs_l42/`).

### 2.2 The new sizes

| m | N | msolve state at the wall | rounds | `D` lower bound | largest matrix | peak RSS | VmPeak |
|--:|--:|:--|--:|--:|:--|--:|--:|
| 3 | 18 | in the 2nd degree-7 step | 23 | 7 | 170,436 × 184,453 | 2.58 GB | 3.66 GB |
| 3 | 21 | in the 2nd degree-7 step | 22 | 7 | 170,445 × 194,842 | 2.58 GB | 3.66 GB |
| 3 | 24 | see `results.csv` (filled in §2.4) | | | | | |
| 4 | 16, 20 | see `results.csv` (filled in §2.4) | | | | | |

- At `N = 18` the committed `f4_fp_tower` run needed six degree-7 steps and
  176 minutes on four unloaded threads. Under this host's load, msolve
  finished one degree-7 step in 45 minutes.
- The first degree-7 step took 2,024 s wall and 384 s CPU at `N = 18`, and
  2,041 s wall and 396 s CPU at `N = 21`.

### 2.3 The trace is independent of `t` up to the first degree-7 step

The `(deg, sel, pairs, rows, cols, new, zero)` rows of `logs/trace_*.txt`:

- `N = 18` and `N = 21` agree on rounds 1–19 in every field, `cols`
  included. In round 20, the first degree-7 step, they agree in every field but
  `cols` (101,390 against 101,435), and both add 2,354 elements with 18,321
  zero reductions. They diverge after that.
- `N = 21` and `N = 24` agree on all 21 rounds **in every field, `cols`
  included**.
- `N = 15` agrees with them in every field through round 7. From round 8 on
  its `sel`, `pairs` and `new` stay within a few units, while `rows` and
  `cols` are smaller: 7% fewer rows and 19% fewer columns in the first
  degree-7 step. That step adds 2,352 elements against 2,354.

### 2.4 XL degree, `δ`, and the Hilbert series

| m | N | `D_XL` (method) | `D_F4` | gap |
|--:|--:|:--|--:|--:|
| 2 | 12 | 16 (NTT `δ = 12`; Macaulay rank: 16) | 5 | 11 |
| 2 | 14 | 18 | 5 | 13 |
| 2 | 16 | 20 | 5 | 15 |
| 2 | 18 | 22 | 5 | 17 |
| 2 | 20 | 24 | 6 | 18 |
| 2 | 22 | 26 | 6 | 20 |
| 2 | 24 | 28 | — | — |
| 2 | 26 | 30 | — | — |
| 3 | 9 | 13 (Macaulay rank, both targets) | 6 | 7 |
| 3 | 12 | ≥ 14 (Macaulay rank; `1 ∉` span at `D = 13`, 40-min cap) | 6 | ≥ 8 |
| 4 | 8 | 12 (Macaulay rank, both targets) | 6 | 6 |

- **`m = 2`.** On every target `S̄` has no zero on the grid (`σ = 0`). The
  coefficients of `1/S̄`, counted by weight, are exactly `C(N, k)` for every
  `k`, so the inverse is fully dense and `δ = N`.
- **Proposition 2 is an equality here.** `NF(I_(D)) = S̄·R_{≤D−4}`, so
  `D_XL = 4 + δ = N + 4` exactly. The direct Macaulay rank at `N = 12`
  reproduces it independently: 1 enters at `D = 16`.
- **The gap grows by one per unit of `N`** at `m = 2`, from 11 to 20. This
  supports H-XL as the protocol fixed it.
- **The affine Hilbert function of the homogenised generators** at `m = 2`,
  `N = 12` is `#sqfree(≤D) − #sqfree(≤D−4)` up to `D = 15`, then 0 at 16
  (`logs/macaulay_*m2_t6*.json`). At `m = 3`, `N = 9` it is still 577 when 1
  enters at `D = 13`.
- **The top-part (degree-of-regularity) series** (`top`):
  - `m = 2`: `top(S) ∈ (y²)`, so the series is `(1+t)^12` and `d_reg = 13`
    at `N = 12`;
  - `m = 3`: the surviving top form is `c·y0_0·y1_0·u0²`. The Hilbert function
    stabilises at 384, the ideal has dimension 1 at infinity, and there is no
    `d_reg`;
  - `m = 4`: the top ideal is not monomial (it has `u0²u1² + …`) and is also
    positive-dimensional.

  No Hilbert-series quantity reproduces `D_F4` at any `N ≤ 18`. By the
  protocol's rule, none is trusted beyond.

**How this resolves the conflict.**
- The `d + δ` bound concerns the Macaulay space `I_(D)`: combinations
  `Σ h_k f_k` with every `deg(h_k f_k) ≤ D`.
- F4 at step degree `D` builds rows `m·g` with `deg(m·g) ≤ D`, where `g` can
  be a *previously reduced* element whose representation in the `f_k` has a
  much higher degree (a degree fall).
- For Kummer towers the falls are everywhere: `y_j·y_j = y_{j+1}` is a carry,
  and §2.3 shows F4 re-deriving the same low window of the tower at every
  `t`. F4's degree is therefore bounded by the *last fall degree*, not by
  `d + δ`.
- The two numbers measure different algorithms. Their gap is now measured
  (11 → 20 at `m = 2`), and the conflict is not evidence against the plateau.
- It is not evidence for the plateau either. A fall-based argument still has
  to show that the degree-7 phase closes without degree 8 as `t` grows, and
  that is what §2.2's wall left open.

## 3. Where the wall falls

- **In this item (5 GB cap, 45 min per job, shared host at load 22–25 on
  4 cores):**
  - msolve finishes `m = 3` at `N ≤ 12` and `m = 4` at `N ≤ 12`.
  - At `m = 3`, `N = 15` it is stopped by **time** inside its second
    degree-7 step, at 1.47 GB.
  - At `N = 18`, 21 and 24 it is stopped by time inside the second degree-7
    step, at 2.6 GB resident.
  - No run reached the memory cap.
- **The memory wall itself was not located.** The second degree-7 step at
  `N = 21` has 170,445 × 194,842 rows and columns. At `N = 18`, `f4_fp_tower`
  needed 13 GB in its degree-7 steps 22–26.
  - Whether msolve fits `N = 21` in 5 GB is unknown: its resident memory was
    still rising when the time cap fired.
  - That is a statement about what was seen, not an extrapolation.

## 4. Deviations from the protocol

1. **The prime.** `p₁ ≈ 2^31` is used, not the design text's "`p ≈ 2^24`",
   because the plateau was measured at `p₁`. This was declared in
   `PROTOCOL.md`.
2. **Controls dropped.** Two control replicates were started and then aborted
   by the operator to free the job slots: `m = 2`, `N = 20` target 1 and
   `m = 3`, `N = 15` target 1. They are kept as `logs/aborted_*.log` and
   `progress.txt` records them; neither is used. The `m = 2`, `N = 22`
   controls were never started.
3. **Concurrency.** Up to three msolve jobs of this item ran at once. Each was
   under its own 5 GB cap, and their summed resident memory stayed below
   5.5 GB.
4. **Arm PL (planted systems) was not run.** The instrument check rests on 15
   agreeing verdicts and degrees instead.
5. **Exhaustive verification.** None was needed for new verdicts, because no
   new system finished.

## 5. Reproduce

```bash
cargo build --release --example pkm_tower_pilot --example ic_next10_03_pkm_degree
cd research/ic_next10_20261014/03-pkm-degree
B=../../../target/release/examples/pkm_tower_pilot
X=../../../target/release/examples/ic_next10_03_pkm_degree
# draw the new sizes (stopped at once) and export them
$B --engine tower --p 2013265921 --kinds kummer --controls tower --planted 0 --random 2 \
   --ladder-t none --budget 60 --max-nnz 1000 --m 3 --t-min 7 --max-t-m3 7 --out gen/m3_t7.jsonl
$B --export-magma-rows gen/m3_t7.jsonl --export-dir gen/magma
$X convert gen/magma/pkm_kummer_p2013265921_m3_t7_x221026520.magma systems/pkm_kummer_p2013265921_m3_t7_x221026520.ms
./run.sh pkm_kummer_p2013265921_m3_t7_x221026520       # msolve -t 1 -v 2, 5 GB, 45 min
./parse.sh > results.csv
$X delta gen/magma/pkm_kummer_p2013265921_m2_t12_x1002906059.magma    # D_XL = N + 4
$X macaulay ../../magma_f4_comparison_20261013/scripts/pkm_kummer_p2013265921_m3_t3_x388393539.magma 40
$X top ../../magma_f4_comparison_20261013/scripts/pkm_kummer_p2013265921_m3_t3_x388393539.magma
python3 -I plot.py
```

## 6. Files

- `PROTOCOL.md`: the protocol, frozen before any data.
- `results.csv`: one row per msolve run.
- `xl.csv`: `D_XL` per system.
- `f4_committed.csv`: the committed `f4_fp_tower` degrees.
- `run.sh`, `parse.sh`, `plot.py`.
- `progress.txt`: the run log, with exits; `124` is the 45-minute cap.
- `runs/`: msolve's `-v 2` logs, outputs and memory samples.
- `runs_l42/`: the linear-algebra-variant check.
- `logs/`: traces, `δ` (`delta_m2.jsonl`), Macaulay and top-part output, and
  the conversion index.
- `systems/`: msolve inputs.
- `gen/`: the stopped draws and their Magma exports.
- `degree_vs_N.svg` and `degree_vs_N.png`: the graph.
- `README.pdf`: this report.

![Solving degree against N](degree_vs_N.png)

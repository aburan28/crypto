# Lead 16.1: where kernel saturation pays, the feasibility crossover in m

Round 4 of `RESEARCH_ISOGENY_FIBERS_ECDLP.md` (§16.1). Toy Koblitz curves
only. Nothing here bears on a deployed curve, and nothing here is a sub-rho
claim.

## Question, hypothesis, falsifier

**Question.** The affine law (Prop. 15.2) charges a base closed under the
rational 2-torsion kernel `M = ⟨T⟩` a mitm cost factor `|M|^{m−2}`. The
u-frame (`u = 1/(x+1)`, `F_u = {u ∈ V}`, closed under `+T` because `1 ∈ V`)
buys symmetrised summation polynomials of lower degree (S₄: boolean degree 4
in `m(ℓ−1)+1` unknowns, against degree 6 in `mℓ` for the direct x-frame S₄).
At what `(n, m)` is the plain x-frame arm infeasible under a shared cap
while the symmetrised arm is not? And does such an arm then beat the best
feasible arm (m = 2, mitm, rho) on the same curve?

**Hypothesis (law, §15.4 algebraic-oracle row).** A crossover cell may exist,
but its cost per relation and total cost stay above the best feasible arm at
the same `n`.

**Falsifier (stated before the runs).** A crossover arm whose total cost
beats every feasible arm at the same `n` would be a counterexample to the
§15.4 row.

## What ran

All computation is Rust. Python only plots the recorded numbers.

**A. Per-attempt panel.** `examples/r4_16_1_m_crossover.rs`. It calls the
library's oracle functions with engine `InheritedF4 { max_degree: 3 }`, the
default split rule, and **one counted cap: `node_budget = 4096` reductions
per solve** for every arm. The driver refuses to start if any engine
environment knob is set. Each `(n, m)` gets one random non-invariant
`V ∋ 1` of dimension `ℓ = ⌈(n + log₂ m!)/m⌉` (the repository's E1
convention), shared by both frames. Two extra cells use the
Frobenius-invariant `V` of divisor `0;1` (`ℓ = 9` at n = 17, `ℓ = 12` at
n = 23), which is the geometry of the earlier n = 17 finding and of the `ic`
anchors. Each curve has one list of 24 targets, uniform in the prime-order
subgroup, shared by every m and every arm. Every verdict is gated against
exact enumeration on the arm's own base. A found relation must re-sum to
`R`, or to `R` with `T` appended, and a refutation must match enumeration.
The result was **0 gate failures in all 25 cells.**

| m | plain x-frame arms (`F_x = {x ∈ V}`) | kernel-saturated u-frame arm (`F_u = {u ∈ V}`) |
|---|---|---|
| 2 | `x-direct`: S₃ | `sym`: symmetrised S₃ |
| 3 | `x-direct`: S₄, degree 6 (the `ic` engine's `descent-algebraic`); `x-chained`: two S₃ links (`groebner_decompose`) | `sym`: symmetrised S₄, degree 4 |
| 4 | `x-chained`: three S₃ links, `4ℓ+2n` unknowns | `sym-chained`: three symmetrised S₃ links, `4(ℓ−1)+1+2n` unknowns |

There is no direct S₅ in either frame in the library, so m = 4 runs only the
chained systems.

```
bash research/isogeny_fibers_20261010/round4/16_1_m_crossover/run.sh panel
# = r4_16_1_m_crossover --curves 17:1 --ms 2,3,4 --v-spec 17:3:div=0.1 --targets 24 --node-budget 4096 --cell-wall-s 300
#   r4_16_1_m_crossover --curves 19:0 --ms 2,3,4 --targets 24 --node-budget 4096 --cell-wall-s 300
#   r4_16_1_m_crossover --curves 23:0 --ms 2,3,4 --v-spec 23:3:div=0.1 --skip x-direct@3 --targets 24 --node-budget 4096 --cell-wall-s 300
```

**B. End-to-end anchors.** These go through `ic bench`, the repository's
index-calculus engine. Each run recovers a planted logarithm, verifies it,
and runs 8 matched counted rho runs. Seed 123212651130. The bases are the
invariant `koblitz-orbit` / `koblitz-symmetrised` bases on divisor `0;1`;
n = 19 has no invariant `V ∋ 1` other than `F₂`, so it uses
`binary-subspace:dimension=10`. The node budget is the same 4096. The m = 3
arms are capped at 300 trials and 10 minutes of contended wall time.
`run.sh anchors` holds the exact commands; each `results/*.log` repeats its
own command line.

Curves: `icv1-f2m17-tm101-00378d4e` (a = 1, r = 65 587, h = 2),
`icv1-f2m19-t797-b6cf2467` (a = 0, r = 130 873, h = 4),
`icv1-f2m23-t5197-69e76b73` (a = 0, r = 2 095 853, h = 4). Host details are
in `results/host.txt`: a 4-core Xeon shared with about 6 other agents, so
every wall time is **contended**.

**Unit.** The panel counts **word XORs** charged by the engine
(`f4_word_ops_thread` deltas), which cover elimination, specialisation and
linear elimination. The Macaulay build and the symbolic substitutions are not
counted, so every oracle cost is a lower bound. `S_rel` is derived as
XOR/relation × (columns + 1) × (GAE per XOR) / √r. The GAE-per-XOR factor is
the measured ns/XOR (`calibrate_word_xor`) over the measured ns per
single-word addition: 0.0035, 0.0049 and 0.0032 at n = 17, 19, 23. `S_rel` is
a **stage diagnostic** for the relation phase only. The end-to-end S is in
table B.

**Feasible** (rule fixed in the driver before the runs): at least one
relation found, (found + refuted)/attempts ≥ 0.5, and zero gate failures.

## A. Measured panel (24 targets per cell, node budget 4096)

`hits` is the number of the 24 targets that exact enumeration says decompose
on that base (measured). For `sym` the count includes `R + T`.

| n | m | ℓ | V | arm | vars | \|F\| | cols | found / refuted / budget | hits | XOR per relation | S_rel (derived) | feasible | wall s (contended) |
|--:|--:|--:|:--|:--|--:|--:|--:|:--|--:|--:|--:|:--|--:|
| 17 | 2 | 9 | rand | x-direct | 18 | 483 | 241 | 17/7/0 | 17 | 4.10e5 | 1.37e3 | yes | 0.07 |
| 17 | 2 | 9 | rand | sym | 17 | 497 | 124 | 9/15/0 | 9 | 6.56e5 | 1.14e3 | yes | 0.11 |
| 17 | 3 | 7 | rand | x-direct | 21 | 129 | 63 | 12/0/12 | 23 | 5.49e8 | 4.87e5 | yes (0.50) | 26.5 |
| 17 | 3 | 7 | rand | x-chained | 38 | 129 | 63 | 23/1/0 | 23 | 2.42e6 | 2.15e3 | yes | 0.65 |
| 17 | 3 | 7 | rand | sym | 19 | 149 | 37 | 12/0/12 | 18 | 5.77e6 | 3.04e3 | yes (0.50) | 4.3 |
| 17 | 3 | 9 | **inv** | x-direct | 27 | 443 | 221 | 0/0/1 (wall guard) | 1 of 1 | — (1.25e9 XORs on one target) | — | **no** | 905.6 |
| 17 | 3 | 9 | **inv** | x-chained | 44 | 443 | 221 | 24/0/0 | 24 | **1.33e6** | 4.10e3 | yes | 0.26 |
| 17 | 3 | 9 | **inv** | sym | 25 | 409 | 102 | 14/0/10 | 24 | **2.26e7** | 3.23e4 | **yes (0.58)** | 5.3 |
| 17 | 4 | 6 | rand | x-chained | 58 | 53 | 26 | 22/0/2 | 22 | 2.36e7 | 8.82e3 | yes | 7.8 |
| 17 | 4 | 6 | rand | sym-chained | 55 | 53 | 13 | 1/0/23 | 5 | 2.75e8 | 5.33e4 | no | 4.6 |
| 19 | 2 | 10 | rand | x-direct | 20 | 1079 | 537 | 13/11/0 | 13 | 1.39e6 | 1.01e4 | yes | 0.18 |
| 19 | 2 | 10 | rand | sym | 19 | 985 | 245 | 10/14/0 | 10 | 1.50e6 | 4.97e3 | yes | 0.17 |
| 19 | 3 | 8 | rand | x-direct | 24 | 257 | 127 | 9/0/15 | 24 | 7.12e9 | 1.23e7 | no | 92.5 |
| 19 | 3 | 8 | rand | x-chained | 43 | 257 | 127 | 24/0/0 | 24 | 4.65e6 | 8.00e3 | yes | 1.0 |
| 19 | 3 | 8 | rand | sym | 22 | 245 | 61 | 6/0/18 | 18 | 3.45e7 | 2.88e4 | no | 7.3 |
| 19 | 4 | 6 | rand | x-chained | 62 | 87 | 42 | 23/0/1 | 24 | 2.72e7 | 1.57e4 | yes | 7.3 |
| 19 | 4 | 6 | rand | sym-chained | 59 | 49 | 12 | 0/0/24 | 1 | — | — | no | 5.2 |
| 23 | 2 | 12 | rand | x-direct | 24 | 4055 | 2024 | 15/9/0 | 15 | 5.86e6 | 2.58e4 | yes | 0.63 |
| 23 | 2 | 12 | rand | sym | 23 | 4209 | 1051 | 3/21/0 | 3 | 3.96e7 | 9.07e4 | yes | 0.93 |
| 23 | 3 | 9 | rand | x-direct | 27 | 467 | 232 | not run (`--skip`, see Limits) | 20 | — | — | — | — |
| 23 | 3 | 9 | rand | x-chained | 50 | 467 | 232 | 9/0/15 | 20 | 7.53e7 | 3.82e4 | no | 5.6 |
| 23 | 3 | 9 | rand | sym | 25 | 497 | 124 | 1/0/23 | 11 | 5.06e8 | 1.38e5 | no | 9.8 |
| 23 | 3 | 12 | **inv** | x-direct | 36 | 4235 | 2116 | not run (`--skip`); `ic` anchor: no output in 10 min | 24 | — | — | — | — |
| 23 | 3 | 12 | **inv** | x-chained | 59 | 4235 | 2116 | 23/0/1 | 24 | **1.03e7** | 4.76e4 | yes | 1.3 |
| 23 | 3 | 12 | **inv** | sym | 34 | 4049 | 1012 | 1/0/23 | 24 | 4.57e9 | 1.01e7 | no | 23.0 |
| 23 | 4 | 7 | rand | x-chained | **74** | 127 | 62 | engine layout refuses (> 64 unknowns) | 17 | — | — | no | — |
| 23 | 4 | 7 | rand | sym-chained | **71** | 137 | 34 | engine layout refuses (> 64 unknowns) | 5 | — | — | no | — |

Columns are the distinct `±[h]P` classes. On `F_u`, P, −P, P+T and −P+T
share a class (Prop. 4.1), so `cols_u ≈ |F_u|/4`. On a random `V` there is no
Frobenius fold. The invariant cells' columns are counted the same way, before
the n-fold that the `ic` bases apply.

**Yield, measured against the law.** On the same `V`, exact enumeration
finds these decomposable-target counts, u-frame against x-frame:

| m | n = 17 | n = 19 | n = 23 |
|--:|:--|:--|:--|
| 2 | 9 vs 17 | 10 vs 13 | 3 vs 15 |
| 3 | 18 vs 23 | 18 vs 24 | 11 vs 20 |
| 4 | 5 vs 22 | 1 vs 24 | 5 vs 17 |

The m = 4 loss is the largest, as `|M|^{m−1}` predicts. At m = 3 the
u-frame's 24/24 on the invariant `ℓ = 9` and `ℓ = 12` bases saturates the
count and carries no information.

## B. Measured end-to-end anchors (`ic bench`, invariant bases)

S = total GAE / √r over the whole pipeline. The table gives the mean over the
repeats, with the range in brackets. Every completed run recovered and
verified the planted logarithm.

| curve | arm | repeats | cols | hit rate | trials / rels | solver budget hits | S | S / S_rho |
|:--|:--|--:|--:|--:|:--|--:|--:|--:|
| icv1-f2m17-tm101-00378d4e | rho (negation map, 8 runs) | 8 | — | — | — | — | 1.76 | 1 |
| icv1-f2m17-tm101-00378d4e | mitm-frobenius, m = 2 | 3 | 14 | 0.60 | 53 / 31 | — | **17.1** (17.0–17.2) | 9.7× |
| icv1-f2m17-tm101-00378d4e | x descent-algebraic, m = 2 | 3 | 14 | 0.60 | 57 / 33 | 0 | 757 (492–999) | 429× |
| icv1-f2m17-tm101-00378d4e | sym, m = 2 | 3 | 13 | 0.32 | 67 / 20 | 0 | 3 112 (2 423–3 798) | 1 764× |
| icv1-f2m17-tm101-00378d4e | **sym, m = 3** | 1 | 13 | 0.50 | 10 / 5 | 5 | **74 460** | 46 567× (vs its matched rho 1.60) |
| icv1-f2m17-tm101-00378d4e | **x descent-algebraic, m = 3** | 1 | 14 | — | — | — | **no output: wall cap 10 min (exit 124)** | — |
| icv1-f2m19-t797-b6cf2467 | rho | 8 | — | — | — | — | 1.63 | 1 |
| icv1-f2m19-t797-b6cf2467 | mitm, m = 2 (binary-subspace, dim 10) | 3 | 502 | 0.59 | 785 / 469 | — | **1 398** | 860× |
| icv1-f2m19-t797-b6cf2467 | x descent-algebraic, m = 2 (dim 10) | 3 | 502 | 0.59 | 785 / 469 | 0 | 34 766 (19 831–46 251) | 21 390× |
| icv1-f2m23-t5197-69e76b73 | rho | 8 | — | — | — | — | 1.09 | 1 |
| icv1-f2m23-t5197-69e76b73 | mitm-frobenius, m = 2 | 3 | 94 | 0.63 | 201 / 126 | — | **143** | 131× |
| icv1-f2m23-t5197-69e76b73 | x descent-algebraic, m = 2 | 3 | 94 | 0.62 | 282 / 175 | 0 | 16 815 (16 218–17 135) | 15 358× |
| icv1-f2m23-t5197-69e76b73 | sym, m = 2 | 3 | 89 | 0.35 | 306 / 108 | 0 | 13 830 (12 068–15 665) | 12 631× |
| icv1-f2m23-t5197-69e76b73 | **sym, m = 3** | 1 | 89 | 0.22 | 266 / 59 | **207** | **938 966** | 921 513× (vs matched 1.02) |
| icv1-f2m23-t5197-69e76b73 | **x descent-algebraic, m = 3** | 1 | 94 | — | — | — | **no output: wall cap 10 min (exit 124)** | — |

The `ic` engine has no chained-x oracle for Koblitz bases. Its x arm at
m = 3 is the direct S₄ descent only.

## Verdict: **confirmed (the law's prediction holds); the falsifier is not met**

1. **A crossover cell exists only against the direct S₄ descent.** The
   smallest such cell is **(n = 17, m = 3)** on the invariant `ℓ = 9` base.
   There the direct x-frame arm is infeasible (panel: one target ran 905 s to
   its node budget; `ic`: no output in 10 min), while the symmetrised arm is
   feasible (panel: 14/24 found, 10 budget, 2.26e7 XOR per relation; `ic`:
   logarithm recovered, S = 74 460). This reproduces the earlier n = 17
   observation. The same holds at n = 23 on the `ℓ = 12` base, where `ic` sym
   m = 3 recovered the logarithm at S = 938 966 and the direct arm returned
   nothing in 10 min. In the panel the sym arm on that base is budget-limited
   on 23 of 24 targets, so by the panel's rule it is infeasible too. On the
   E1-sized random `V` there is no strict crossover: at n = 17 both arms sit
   at exactly 0.5 solved, and at n = 19 both are infeasible.
2. **Against the best plain arm there is no crossover anywhere in
   n ∈ {17, 19, 23}, m ∈ {3, 4}.** The chained x-frame system, which is
   already in the library (`groebner_decompose`), is feasible wherever the
   symmetrised arm is. It is cheaper per relation by 2.4× (n = 17, ℓ = 7),
   17× (n = 17, ℓ = 9 inv), 7.4× (n = 19, ℓ = 8) and 440× (n = 23, ℓ = 12
   inv). At m = 4 the kernel-saturated chain is the infeasible one
   (1/24 and 0/24 found, the rest budget-limited), and the x chain solves
   22/24 and 23/24.
3. **The crossover arm loses to every feasible m = 2 arm and to rho.** At
   n = 17, `ic` gives sym m = 3 S = 74 460, against x m = 2 S = 757 (98× less),
   mitm m = 2 S = 17.1 (4 360× less) and rho 1.76. At n = 23 the gaps grow: sym
   m = 3 S = 938 966 against mitm m = 2 S = 143 (6 570×). The relation-phase
   diagnostic agrees: at n = 17 the sym m = 3 cell has S_rel 3.2e4, against
   1.4e3 for x-direct m = 2 and 4.1e3 for x-chained m = 3 on the same base.

So kernel saturation does buy feasibility, but only against the single plain
system the `ic` engine happens to plug in (direct S₄). It never beats the
best feasible arm at the same n, which is what §15.4's algebraic-oracle row
predicts.

## Limits (what was not done)

- **x-direct at n = 23, m = 3 was not run in the panel.** In the first
  attempt the process ran past the 30-minute background limit while still
  inside the first `ℓ = 9` target, and was killed
  (`results/panel_n23_b4096_killed_attempt1.log`). The rerun skips that arm
  (`--skip x-direct@3`) and relies on the `ic` anchor, which produced no
  output in 10 min, and on n = 17, where one 27-unknown solve ran 905 s.
  "Infeasible" there is a wall-time observation under the shared counted cap,
  not a counted result. The per-cell wall guard checks only between targets,
  so a single solve can run past it.
- **No direct S₅** in either frame. The library has none in characteristic 2,
  so m = 4 compares chained systems only, and the symmetrised chain is not the
  [2,2,2,2,1]-type degree reduction the lead names. A direct symmetrised S₅
  (182 unknown F₂ coefficients, A + s·B with s² = Σw + s) is the natural next
  arm. I did not build it.
- **m = 4 at n = 23 is out of reach for both arms.** The chained layouts need
  74 and 71 unknowns, and the engine caps at 64 (`MAX_VARS`).
- 24 targets per cell. The hit rates carry binomial noise of roughly ±0.1.
  The `ic` m = 3 rows are single repeats (10 and 266 trials).
- **Panel against `ic`, sym m = 3 on the n = 23 `ℓ = 12` base.** The panel
  found 1/24 with 23 budget-limited. `ic` found 59 relations in 266 trials,
  with 207 budget-exceeded calls. The rates differ (0.04 against 0.22). I did
  not resolve this; the target sampling is a candidate cause. Both runs say
  the arm is mostly budget-limited.
- The XOR counts exclude the Macaulay build and the substitutions. The
  word-XOR-to-GAE price was measured once per process on a contended host.
- The panel's S_rel prices only the relation phase, with cols + 1 relations.
  It is not an end-to-end S. Table B is.
- Conductor was not installed in this environment, so `conductor check` could
  not run.
- `cargo clippy --release --example r4_16_1_m_crossover -- -D warnings` fails
  on warnings that already exist in `src/` (29 library warnings, including an
  unknown-lint E0602). The example adds no warnings of its own: its single
  warning is that duplicate. `rustfmt --check` passes.

## Next step

Build the direct symmetrised S₅ (`A(w) + s·B(w)`, symmetric in five `w`'s)
and the direct x-frame S₅ by interpolation. Run both at n = 17 and 19 with
`ℓ = 6` under the same budget. That cell is where the lead's degree
reduction actually lives at m = 4. The prediction from this round is that
the x chain still beats it.

## Files

- `examples/r4_16_1_m_crossover.rs`: the panel driver.
- `run.sh`: exact commands for A and B.
- `results/panel_n{17,19,23}_b4096.json` and `.log`: the panel, including the
  per-attempt records. `results/panel_n23_b4096_killed_attempt1.log` is the
  killed first attempt.
- `results/{k1_17,k0_19,k0_23}__*.json` and `.log`: the `ic bench` anchors.
  The x-m3 runs have only a `.log`, with exit 124.
- `results/host.txt`.
- `plot.py`, `crossover_data.csv`, `crossover.svg`, `crossover.png`: the
  graph.
- `report.typ`, `report.pdf`.

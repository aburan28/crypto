# A supersingular auxiliary curve for PKM's isogeny tower (2026-10-15)

Does PKM's isogeny-tower factor base need the ordinary auxiliary curve whose
P-192 search the backlog card
[`XP-20261015-p192-isogeny-aux-search-t48`](../../backlog/XP-20261015-p192-isogeny-aux-search-t48/README.md)
prices at about 11.2 core-years, or does the supersingular curve
`y² = x³ + x` (one draw at P-192) do the same job? Three pre-registered
questions: the structure of `V` (Q1), the F4 solving degree (Q2) and the
charged relation phase (Q3).

- **Protocol.** [`PROTOCOL.md`](PROTOCOL.md), committed in `8d719fddf` before
  any measured cell (sha256 `db109df9…c1ad65`). Its disclosure section lists the
  scratch dry runs made before it.
- **Code.** `src/cryptanalysis/ic_framework/prime_tower.rs`
  (`supersingular_tower`, `prime-isogeny:…,aux=supersingular`),
  `examples/pkm_ss_tower.rs` (`structure`, `f4`),
  `examples/prime_tower_curves.rs --plus`.
- **Runs.** [`run_q1.sh`](run_q1.sh), [`run_q2.sh`](run_q2.sh),
  [`run_q3.sh`](run_q3.sh); [`progress.txt`](progress.txt) records every start,
  end and exit status, and the two container restarts.
- **Tables.** `q1/pooled.tsv`, `q1/vs_torus.tsv`, `q1/positive_control.tsv`
  (`q1/tabulate.sh`); `q2/cells.tsv` (`q2/tabulate.sh`);
  `q3/tables/summary-arms.tsv`, `q3/tables/summary-ratios.tsv`
  (`q3/tabulate.sh`; the ratios and intervals are ecbench's).
- **Graphs.** [`figures/plot.py`](figures/plot.py) draws
  `figures/q1_equivalence.svg`, `figures/q2_d_vs_n.svg` and
  `figures/q3_s_and_ratios.svg` (PNG beside each) from those TSVs.
- **Report.** [`report.typ`](report.typ) → [`report.pdf`](report.pdf); run
  report [`reports/run-20261015T2315Z/`](reports/run-20261015T2315Z/report.pdf).

**Labels.** *exact* = a deterministic, verified computation; *measured* = a
value in a committed row or sealed ecbench record; *derived* = an argument
or a calculation shown here.

## 0. Verdicts

| question | verdict | class |
|:--|:--|:--|
| **Q1**: is `V_ss` the Dickson torus set up to a Möbius change of `x`? | **No; by the letter of the rule *partial*.** Over 103 primes `p ≡ 3 mod 4` (103–2999) and `t = 3, 4, 5`, none of the 6,507 `y² = x³ + x` sets and none of the 6,768 `j = 0` sets (both twists) is Möbius-equivalent to any torus set; 2 of the 12,969 `y² = x³ − x` sets are, at one cell (`p = 383, t = 3`, one shared class), which the pre-registered rule counts as *partial* rather than H1_distinct. Ordinary curves: 0 of 24,580. At `t = 2` every family, ordinary included, matches 75 % of the time (one modulus). So the supersingular route does not inherit round 8's verdict through Q1. | exact |
| **Q2**: F4 solving degree, supersingular vs ordinary vs torus | **H_same fails; H2a holds, H2b fails.** The supersingular tower gives the ordinary isogeny tower's `D` in every one of the 14 `(m, N)` cells and at all three primes, for all four supersingular curves (`y² = x³ ± x`, `y² = x³ + 1`, `y² = x³ + d³`): `D = 3, 5, 5, 5, 5, 5, 5, 5, 5` at `m = 2, N = 2…18`, `6` at `N = 20`, with widths equal to the ordinary tower's to within one column. The Dickson torus gives `4` at `N = 4…10` (and `5` vs `6` at `m = 3, N = 6`). | measured |
| **Q3**: relation phase, charged | **The search cost disappears and nothing is lost per base element.** With `aux=supersingular` the auxiliary work is 0.9–1.9 % of `S` (the build alone) against 10–27 % for the searched ordinary curve; `S` is within 1.3–2.2 % of a plain base of the same size, and the relation yield is unchanged. On s22 and s24 the supersingular arm costs 0.82× and 0.74× the ordinary arm; on s20 it costs 1.11×, because its base is 15 % larger (141 against 123 abscissae; the plain bases of the two sizes differ by 1.21×). H3a, H3d and both controls hold; H3b and H3c fail on s20 by the letter of the rule. | measured |

**What this means for the P-192 card (derived).** At these sizes a
supersingular `E'` is indistinguishable from an ordinary one in the solving
degree of the decomposition system (Q2), and it removes the search from the
relation phase without changing the per-base-element cost (Q3). The one
structural difference found is that neither family is the Dickson torus
family (Q1, Q2), and the supersingular route shares that with the ordinary
one. So the 11-core-year ordinary search buys nothing these experiments can
see: a supersingular `E'` (one draw) is a sufficient input for the relation
phase, unless a reason outside this study requires an ordinary curve. This
does not touch the card's real gate, lead 1 (whether `D` stays bounded as `N`
grows), and says nothing about P-192's security.

## 1. Q1 — the structure of `V` (exact)

`pkm_ss_tower structure`, primes `p ≡ 3 mod 4` with `100 ≤ p ≤ 3000`,
`t = 2…5` wherever `t + 1 ≤ v₂(p + 1)`. For each `(p, t)` it enumerates every
tower set of every family and compares their complete `PGL₂(F_p)` invariants
(the least image of the set under the Möbius maps sending an ordered triple
of its points to `(0, 1, ∞)`).

| `t` | cells `(p, t)` | `y² = x³ + x` → torus | `y² = x³ − x` → torus | `y² = x³ + 1` → torus | `y² = x³ + d³` → torus | ordinary (3 curves) → torus | positive control → `y² = x³ + x` |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 2 | 103 | 13,965 / 18,617 | 27,883 / 37,181 | 7,221 / 9,638 | 7,186 / 9,638 | 54,738 / 72,689 | 18,617 / 18,617 |
| 3 | 53 | 0 / 4,993 | **2** / 9,959 | 0 / 2,593 | 0 / 2,593 | 0 / 18,763 | 4,993 / 4,993 |
| 4 | 27 | 0 / 1,256 | 0 / 2,500 | 0 / 650 | 0 / 650 | 0 / 4,944 | 1,256 / 1,256 |
| 5 | 12 | 0 / 258 | 0 / 510 | 0 / 141 | 0 / 141 | 0 / 873 | 258 / 258 |

Entries: sets equivalent to some set of the target family / all sets, pooled
over the cells (`q1/pooled.tsv`; per cell in `q1/vs_torus.tsv`). The `j = 0`
columns cover the 52 / 26 / 13 / 6 cells with `p ≡ 2 mod 3`. The torus family
has 18,617, 4,993, 1,256 and 258 sets at `t = 2…5`.

**The one exception.** At `p = 383`, `t = 3`, one Möbius class is shared by a
torus set and two `y² = x³ − x` sets. No other `(p, t ≥ 3)` cell has any match for any supersingular
family, and `y² = x³ + x`, the P-192 curve, has none at all. The rule of
PROTOCOL.md (H1_distinct needs zero matches) therefore reads *partial*; the
match is isolated, not family-wide, and is not explained here.

**Other equivalences found (expected).** Ordinary control curves with equal
`j` (isomorphic: `p = 127, t = 4`, `j = 1`; `p = 1151, t = 3`, `j = 24`;
`p = 1663, t = 4`, `j = 32`) share classes, as the positive control predicts.
No supersingular family matches another (`y² = x³ + x` against
`y² = x³ − x`: 0 at `t ≥ 3`).

**Checks.** The positive control (`y² = x³ + 16x`, isomorphic to `y² = x³ + x`
by `(x, y) ↦ (4x, 8y)`) matched 100 % of its sets in all 195 `(p, t)` cells: **passed**. The canonical form agreed with a brute-force
search over all of `PGL₂(F_p)` on all 906 checked (set, family) pairs at `p ≤ 131`. The
counts of `y² = x³ + x` cosets equal `(p + 1 − 2^{t+1})/2^{t+1}` (unit test).

**Why (derived).** A Möbius map taking one fibre of a degree-`2^t` map onto a
fibre of another does not, by itself, relate the maps, so the computation is
needed; but there is a structural reason to expect no family-wide map. The
Dickson set is a fibre of `D_{2^t}`, a polynomial, totally ramified over `∞`.
The supersingular set is a fibre of the composed Vélu x-map `Φ` of a cyclic
`2^t`-isogeny, whose fibres are `x(P + K)`; two points of a fibre coincide only
when `2P ∈ K`, so every ramification index of `Φ` is at most 2, and `Φ` is not
Möbius-conjugate to `D_{2^t}` for `t ≥ 2`. At `t = 2` a 4-point set has one
modulus, and every pair of families (ordinary ones included) shares classes
by coincidence; from `t = 3` on the computation finds one isolated
coincidence (above) and nothing family-wide. Supersingularity does not make
`V` a torus set: `y² = x³ + x` and the torus both have `p + 1` rational points,
but one is a genus-1 group and the other a conic, and no algebraic map carries
cosets of one onto cosets of the other. The supersingular route therefore does
**not** inherit round 8's `H_same` verdict through Q1; Q2 measures it
directly.

## 2. Q2 — F4 solving degree (measured)

Engine `f4_fp_tower` (library), reduced presentation, full-rank exit, degree
bound `n + d + 6`; two random targets per cell; one thread; 336 systems, none
timed out or stopped for size; every target was inconsistent (no
decomposition), and the 296 tower rows with `m ≤ 3` agree with an exhaustive
recount over `V^m` (the 12 `m = 4` rows have no recount; the 28 null rows are
not summation systems).

`D / width / steps`, both targets equal in every cell:

### `m = 2`

| `N` | `ss-j1728` (T24, T31, S25) | `ss-j1728-twist` | `ss-j0`, `ss-j0-twist` (S25) | `ordinary` (T24, T31, S25) | `torus` (T24, T31, S25) | null on `ss-j1728` (T31) |
|--:|:--|:--|:--|:--|:--|:--|
| 2 | 3 / 4 / 2 | 3 / 4 / 2 | 3 / 4 / 2 | 3 / 4 / 2 | 3 / 4 / 2 | 3 / 4 / 2 |
| 4 | 5 / 16 / 4 | 5 / 16 / 4 | 5 / 16 / 4 | 5 / 16 / 4 | **4** / 16 / 4 | 5 / 16 / 4 |
| 6 | 5 / 57 / 7 | 5 / 57 / 7 | 5 / 57 / 7 | 5 / 57 / 7 | **4** / 56 / 6 | 5 / 57 / 7 |
| 8 | 5 / 150 / 9 | 5 / 150 / 9 | 5 / 150 / 9 | 5 / 150 / 9 | **4** / 152 / 8 | 5 / 150 / 9 |
| 10 | 5 / 346 / 12 | 5 / 346 / 12 | 5 / 346 / 12 | 5 / 346 / 12 | **4** / 336 / 10 | 5 / 346 / 11 |
| 12 | 5 / 1182 / 16 | 5 / 1181 / 16 | 5 / 1181 / 16 | 5 / 1181 / 16 | 5 / 1223 / 13 | 5 / 1182 / 15 |
| 14 | 5 / 2629 / 19 | 5 / 2629 / 19 | 5 / 2629 / 19 | 5 / 2629 / 19 | 5 / 2303 / 15 | 5 / 2628 / 18 |
| 16 | 5 / 5034 / 22 | 5 / 5034 / 22 | 5 / 5034 / 22 | 5 / 5034 / 22 | 5 / 4786 / 17 | 5 / 5034 / 21 |
| 18 | 5 / 9220 / 24 | 5 / 9221 / 24 | 5 / 9221 / 24 | 5 / 9221 / 24 | 5 / 8639 / 20 | 5 / 9184 / 23 |
| 20 (T31) | 6 / 35157 / 28 | — | — | 6 / 35158 / 28 | 6 / 33426 / 24 | — |

Every entry is the same at every prime where the arm ran. Positive control
at `p₀ = 786433`: `ordinary` 3, 5, 5, 5, 5 against `split` Dickson 3, 4, 4, 4, 4
(`N = 2…10`), as round 8's R0 — **passed**.

### `m = 3, 4` (`chain`, T31)

| `m` | `N` | `ss-j1728` | `ordinary` | `torus` | null on `ss-j1728` |
|--:|--:|:--|:--|:--|:--|
| 3 | 6 | 6 / 168 / 10 | 6 / 200 / 10 | **5** / 164 / 10 | 6 / 168 / 10 |
| 3 | 9 | 6 / 1196 / 18 | 6 / 1196 / 18 | 6 / 1184 / 16 | 6 / 1200 / 18 |
| 3 | 12 | 6 / 5088 / 25 | 6 / 5086 / 25 | 6 / 5018 / 18 | 6 / 5065 / 24 |
| 4 | 8 | 6 / 959 / 20 | 6 / 959 / 20 | 6 / 1080 / 21 | 6 / 950 / 20 |
| 4 | 12 | 7 / 9453 / 34 | 7 / 9454 / 34 | 7 / 9451 / 29 | 7 / 9387 / 32 |

**The rules of PROTOCOL.md, applied.**
- **H2a (supersingular = ordinary): *same*.** Equal `D` on both targets in
  every cell (m = 2: 9 cells at each of three primes plus `N = 20`; m = 3, 4: 5
  cells), for each of the four supersingular curves. Widths agree to within one
  column and the step counts agree; at `N = 20` the multiply-adds are
  437,488,141,950 and 437,488,139,872 (supersingular, two targets) against
  438,029,444,285 (ordinary, both targets).
- **H2b (supersingular = torus): *differs*.** 5 against 4 at `m = 2`,
  `N = 4…10`, on both targets at every prime, and 6 against 5 at `m = 3, N = 6`.
- **H_same: fails** (it needs both). This is the prediction the protocol wrote
  down.
- **Cross-check with round 8.** The torus rows reproduce round 8's
  `D / width / steps` in every shared cell, and at `N = 20` its multiply-add
  count, 478,782,978,864, equals round 8's split-Dickson count to the unit. The
  ordinary isogeny `N = 20` cell (at `T31`) reproduces round 2's I0 (at `p₀`) (`D = 6`, 35,158
  columns, 28 steps).
- **Null (gates nothing).** `S₃`'s support with random coefficients over the
  supersingular tower gives the supersingular `D` in every cell: as in round 8,
  the degree belongs to the tower algebra and the `S₃` support, not to the curve.

**Why (derived).** A supersingular tower is a Vélu tower like the ordinary
one: its equations `y_{j+1}(y_j − ξ_j) = y_j² − ξ_j y_j + τ_j` have the top
form `y_j y_{j+1} − y_j²`, while the Dickson equations have `−y_j²`
(PKM note §3.4). The F4 computation sees the leading forms and the supports,
which are the same for every Vélu tower whatever `ξ_j, τ_j`; the one-column
difference at `N = 12, 18` comes from `ξ₀ = 0` on `y² = x³ + x` (the kernel
point `(0, 0)` drops a monomial).

## 3. Q3 — the relation phase, charged (measured)

Three curves of the P-192 shape (`v₂(p − 1) = 1`, `v₂(p + 1) = 10`),
registered in this change; one ecbench session each; 7 arms × 8 one-target
workloads × (1 warm-up + 2 measured) rounds; **504/504 runs verified, every
replay identical** (`ecbench verify --replay-all`, audit OK on all three).
L0 (`--cpus none --allow-busy`, a host shared with another lane): operation
counts only, no wall time.

| curve | `p` | `t` | ss base (abscissae) | ordinary base, curves tried | aux share of `S`, ss | aux share, ordinary |
|:--|--:|--:|--:|--:|--:|--:|
| s20 `icv1-fp20-t1143-5bef62bd` | 537,599 | 8 | 141 | 123, 71 | 1.9 % | 10.4 % |
| s22 `icv1-fp22-tm1105-62a9d176` | 2,511,871 | 9 | 246 | 235, 684 | 1.4 % | 24.5 % |
| s24 `icv1-fp24-tm4507-ffb4b7f9` | 11,750,399 | 9 | 253 | 251, 1,043 | 0.9 % | 27.1 % |

Mean `S` (16 measured runs per arm):

| curve | rho (neg.) | ss rootfind | plain (ss size) | ordinary rootfind | plain (ord size) | ss mitm |
|:--|--:|--:|--:|--:|--:|--:|
| s20 | 1.448 | 69.72 | 68.24 | 62.94 | 56.43 | 55.51 |
| s22 | 1.087 | 84.03 | 82.95 | 102.03 | 77.16 | 73.85 |
| s24 | 1.068 | 33.75 | 33.09 | 45.87 | 33.04 | 31.43 |

Paired ratios `Σ S_B / Σ S_A` with ecbench's 95 % cluster-bootstrap interval
(16 pairs on 8 workloads):

| curve | ss ÷ plain (ss size) | ordinary ÷ plain (ord size) | **ss ÷ ordinary** | plain(ss) ÷ plain(ord) | rootfind ÷ mitm (pos.) | A/A rho | ss ÷ rho |
|:--|:--|:--|:--|:--|:--|:--|:--|
| s20 | 1.022 [1.010, 1.035] | 1.115 [1.098, 1.132] | 1.108 [1.092, 1.122] | 1.209 [1.198, 1.220] | 1.256 [1.245, 1.265] | 1.000 | 48.1 [40.2, 57.3] |
| s22 | 1.013 [1.006, 1.019] | 1.322 [1.311, 1.332] | 0.824 [0.820, 0.827] | 1.075 [1.071, 1.079] | 1.138 [1.133, 1.143] | 1.000 | 77.3 [63.9, 98.8] |
| s24 | 1.020 [0.992, 1.048] | 1.388 [1.354, 1.421] | 0.736 [0.715, 0.764] | 1.002 [0.997, 1.005] | 1.074 [1.067, 1.081] | 1.000 | 31.6 [25.0, 41.4] |

Relation yield (`relations_found / targets_tried`, pooled over 16 runs):
s20 0.9995 (ss) / 0.9995 (plain, ss size) / 0.9896 (ordinary); s22 0.9997 /
1.0000 / 1.0000; s24 0.8377 / 0.8526 / 0.8337. Every IC arm's mean `S` is 29–94× rho's with
the negation map (the distance remaining; prime-curve tower IC has not
approached rho at any size measured).

**The rules of PROTOCOL.md, applied.**
- **P1: held.** 504/504 verified; 336/336 replays identical; three audits OK.
- **H3a: held.** `aux_search_gae = 0` on every supersingular run; its aux
  share is 1.9, 1.4 and 0.9 %, all of it the build (`⟨Q⟩`, `V`, the chain,
  and drawing `Q`: 564, 1,085 and 1,089 GAE) that the ordinary arm pays too
  (584, 1,293, 1,314 GAE) on top of its search (2,192, 21,586, 41,211 GAE).
- **H3b: held on s22 and s24, failed on s20.** On s20 the interval
  [1.010, 1.035] excludes 1 and reaches past 1.03. The excess, 2.2 %, matches
  the build's 1.9 % share; the same build cost is inside the ordinary arm's
  `S`. The rule did not allow for it.
- **H3c: held on s22 and s24 (0.82×, 0.74×), failed on s20 (1.11×).** On s20
  the supersingular base has 141 abscissae and the ordinary 123, and the plain
  bases of those sizes differ by 1.21×. Corrected for size by the plain arms,
  `(ss ÷ plain-ss) ÷ (ord ÷ plain-ord)` is 0.92, 0.77 and 0.73 (derived from
  the ratios above, not an ecbench interval). Base size is fixed by `t` and the
  draw of `R'` (about half of `V` lifts to the curve), not by the auxiliary
  curve's type.
- **H3d: held.** The supersingular yield is within 1.8 % of the plain base of
  the same size on every curve.
- **Controls: passed.** A/A rho 1.000 on every curve; rootfind ÷ mitm 1.256,
  1.138, 1.074, above 1 with the interval excluding 1, as in PR #1889.

**Scope.** 20–24-bit curves, one seed per arm, `m = 3` table oracles. The
statement is about the cost of building the base and its per-element cost,
not about any IC exponent; every arm is far above rho.

## 4. What is not done, and caveats

- **No P-192-size run** of anything here: the F4 engine needs `p < 2³²`, and
  the plug-in is `u64`. The P-192 statement is the construction fact (one draw;
  Lane G) plus these small-size equivalences.
- **Q2 covers `N ≤ 20` at `m = 2` and `N ≤ 12` at `m = 3, 4`.** Equal `D` up to
  there does not prove equal `D` beyond; the derived reason (same leading
  forms and supports) is the argument that it continues.
- **Q1 is exact on the primes run**, not a proof for all `p`; the
  ramification argument is the proof sketch that the torus and Vélu maps are
  not conjugate. The `p = 383` coincidence is recorded, not explained.
- **Container restarts.** Two: one before the protocol (no data), and one that
  killed `run_q2.sh` during `T31-ordinary-m34` (its 8 partial rows kept in
  `q2/interrupted/`; the run was repeated whole, and `run_q2.sh`'s 180-minute
  deadline restarted with it). `run_q1.sh` was made resumable in chunks before
  Q1 ran; flags and range unchanged.
- **`fb_id` and the explicit default.** `prime-isogeny:t=8,aux=ordinary`
  builds the same points as `prime-isogeny:t=8` but has another `fb_id`
  (the parameters are hashed); the sessions use the bare form.
- **Registry.** `docs/curves/registry.json` and `src/cryptanalysis/curve_aliases.json`
  gained the three curves (targeted update with the build script's own
  harvester, because the checkout is sparse); the registry-derived artifacts
  (`docs/browser/data.json`, covers, traits, Jacobian transfers, leaderboard)
  are **not** regenerated here.
- Conductor's control plane was unreachable; no task card.

## Environment

Claude Code cloud container, 4 vCPU, 15 GiB, shared with another lane.
Binaries built from `8d719fddf` in this lane's own target directory
(`binaries.sha256`): `ecbench` `92e39467…53b821`, `pkm_ss_tower`
`f277e003…925b23`. One thread, `nice -n 10`, 4 GB address-space cap.

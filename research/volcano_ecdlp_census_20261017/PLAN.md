# Volcano ECDLP census: every curve, every level, three methods, transfer up and down

Status: **preregistered plan, frozen before any data** (2026-10-17). Changes
after Phase B starts go to `AMENDMENTS.md` with a timestamp and a reason.
Outcome-driven changes to workloads, seeds, units or accounting are
inadmissible (ecbench-measure §1).

## 0. Question

For whole ℓ-isogeny volcanoes of prime-field curves:

1. **Census.** List every curve (F_p-isomorphism class, by j and model),
   every verified ℓ-isogeny edge and every level, crater to floor. Store
   them durably in git and on Google Drive.
2. **Cost by level.** On every curve, generate ECDLP instances and measure
   Pollard rho, BSGS and index calculus in ecbench's counted unit `S`
   (group operations) and in wall time. Does any method's cost depend on
   the volcano level, or on the curve within a level?
3. **Transfer.** Carry instances **down** (crater → floor), **up**
   (floor → crater) and **sideways** (along the crater) with verified Vélu
   maps. Measure what transport costs, check that `log` is preserved, and
   find out whether "transport, then solve on the destination" ever beats
   "solve where the instance was born".

Prior evidence it builds on, not repeats: aburan28/crypto#2030 (`volcano-rho`:
no level effect on plain rho; j = 0 / 1728 craters gain from automorphism
folding), aburan28/cryptanalysis#943 (`ca volcano`), #2001 (complete-class
enumeration and graph-wide IC survey), and #2065 (an earlier rho/BSGS/IC volcano
campaign lost its checkpoints because they existed only in the container).
**Lesson from #2065: nothing exists until it is in git and on Drive.** Every
phase ends with an upload and a read-back.

Prior hypothesis (from #2030 and theory): rho and BSGS show **no** level
effect except automorphism folding on special-j craters. IC on a prime field
is not expected to be subexponential, and its cost is not expected to depend on
level. A null result is a result. The success criterion (CLAUDE.md) says to
report IC's distance from matched rho per row; it is not pass/fail.

## 1. Instruments (no new solvers)

| job | tool | where |
|:--|:--|:--|
| curve search (deep volcano, given bits) | `volcano-rho --bits B --depth D [--crater …] [--j0]` | branch of #2030 (`aburan28/youthful-ride-6u51pp`), built in a worktree, search mode only |
| walk, levels, Vélu edge checks, transported `(G, Q)` per vertex | `ca volcano --p --a --b --ell --max-vertices … --diagram … ` | aburan28/cryptanalysis `main` (`build/ca`) |
| rho / BSGS / IC measurement | `ecbench run / verify / compare / table` | this repo, `cargo build --release --bin ecbench` |
| diagram | `tools/volcano-diagram/volcano_diagram.py` | #2030 branch (copy into this study if not merged) |
| independent walk cross-check | `volcano-rho --curve P,A,B` vs `ca volcano` | must report the same vertex count per level |

**Transported instances in ecbench, exactly.** Each curve E' enters a spec as
`prime_explicit` with `gx, gy` = φ(G_src), where φ is the walk path from the
source curve and G_src is the source generator. The spec uses
`"target_kind": "planted_shared"`, so `k_i` depends only on
`(r, target_seed, i)`. Then ecbench's i-th target on E' is
`[k_i]φ(G_src) = φ([k_i]G_src) = φ(Q_src,i)`: the transported instance. φ has
degree ℓ^d with ℓ ≠ r, so it is an isomorphism on the order-r subgroup and keeps
`log`. The executor checks this once per source on ≥ 3 curves: apply `ca`'s
edge maps to `[k_i]G_src` and compare.

Arms (same in every session):

| arm | method | role |
|:--|:--|:--|
| `rho-neg` | `rho.negation` | reference (matched rho) |
| `rho-plain` | `rho.plain` | baseline |
| `bsgs` | `bsgs.textbook` | candidate |
| `bsgs-neg` | `bsgs.negation` | candidate |
| `ic` | `ic.pipeline`, `factor_base` and `oracle` frozen in Phase A from the calibration (default `prime-abscissa:size=32`, `mitm:negation_folded=1`) | candidate |
| `rho-neg-aa` | `rho.negation` | A/A control: wall-time noise floor |
| `rho-aut` (special-j craters only) | the automorphism-folded rho used by #2030 / `ca --automorphisms`, if ecbench exposes it; otherwise via `ca` and labelled "other instrument, ops only" | positive control |

## 2. Volcanoes (population)

Chosen in Phase A by the search tool, using a seed that is recorded before the
search runs. The bit size `b` is fixed by the Phase A calibration rule (§4).

| id | ℓ | root | why |
|:--|:-:|:--|:--|
| V1 | 2 | generic j, depth ≥ 4, p ≡ 1 (mod 8) | main population |
| V2 | 3 | generic j, depth ≥ 3 | second ℓ, different shape |
| V3 | 2 | j = 0 crater (D_K = −3) | special-j positive control |
| V4 (stretch) | 2 | quadratic twist of V1 | the other side of the twist; different order |

"All curves" means **every vertex of each volcano**: every F_p-isomorphism
class reached by ℓ-isogenies, at every level. Each vertex also gets **two
isomorphic re-scalings** (`(u⁴a, u⁶b)`) as the negative control: they are the
same curve in another model, so any "curve effect" must exceed their spread. The
full isogeny class (all ℓ at once, #2001's `volcano-cli`) is listed in Phase A
only, as a census; it is **not** measured.

## 3. Instances ("generators from all sides")

For each volcano, sources are chosen before measuring:

| family | source curve(s) | transport | label |
|:--|:--|:--|:--|
| **D** (down) | the crater root | root → every vertex | `src=crater` |
| **U** (up) | 2 floor vertices in different subtrees, picked by seeded index | floor → every vertex (ascends, crosses the crater, descends other subtrees) | `src=floor-a`, `src=floor-b` |
| **S** (sideways) | a second crater vertex (if crater > 1) | horizontal-first | `src=crater-b` |
| **N** (native) | every vertex | none: a fresh generator per curve (seeded hash to the subgroup) | `src=native` |

Per (vertex, source family): `targets_per_curve = 8`, using a different
`target_seed` per family, all recorded. D and U share `r` and `k_i`, so their
costs pair on the same scalar. This is a check for BSGS (identical counts are
expected and act as a sanity check, not a finding), and a joint resampling unit
for rho and IC.

**Transfer cost.** For every (source, destination) pair, record the path length
`d`, direction counts (up / down / horizontal), Vélu evaluations, and wall time
to map `(G, Q)`. Measure it with `ca`'s point maps, adding a counter in
`volcano.c` if it has none. Report "transport + solve" against "solve at
source" for each method, priced in GAE (field-op estimate per Vélu evaluation
stated, 0 / 0.1 / 1 band as for any unpriced work).

## 4. Phases

Each phase is one executor-subagent assignment. The coordinator (the session
that wrote this plan) reviews each phase's gate before it launches the next.

### Phase A: setup, census, calibration (gate: G-A)
1. Build `ca` (cryptanalysis `main`), `ecbench` (this branch), and
   `volcano-rho` (#2030 branch, worktree under the scratchpad). Record the
   commits.
2. **Calibration.** On roots found at b ∈ {18, 20, 22, 24, 26} bits (V1
   shape), run one ecbench session per b with all arms, 8 targets, 3 rounds.
   Record per-solve mean and CV of `S` and of wall time per arm.
   **Rule:** b = the largest size where the IC arm's mean wall time is ≤ 0.5 s
   per solve **and** every IC run verifies. Also fix the IC factor-base size
   and oracle there, as the best verified IC arm in `S` among ≤ 4 configs
   pre-listed in `calibration/configs.json` before running. Freeze both in
   `FROZEN.json`.
3. Search V1–V4 at that b with recorded seeds. Walk each with `ca volcano`
   (`--max-vertices` large enough for the whole volcano; refuse a truncated
   walk). Cross-check vertex counts per level with `volcano-rho --curve`.
4. Write the census: `census/<V>/vertices.jsonl` (one line per vertex: id,
   level, p, a, b, j, #E, r, cofactor, trace, End-order level, parent edge,
   kernel, controls), `edges.jsonl` (each verified edge, direction, kernel
   x-coords), `transport/<src>.jsonl` (φ(G_src) on each vertex, path),
   `diagram.json/.svg/.png`. Add the all-ℓ class census from `volcano-cli` if
   it builds; otherwise record "not built" and why.
5. Generate every Phase B spec (`specs/<V>-<src>-<chunk>.json`) and the run
   manifest `runs/manifest.json`. Size them with §5 so each session is ≤ 30 min.
   `ecbench plan` must succeed on each spec.
6. **Upload**: commit + push; tar.zst the census and specs, upload to the Drive
   folder, read back the size and the SHA-256 of the first and last MiB, and
   record the file IDs in `DRIVE.md`.

**G-A:** builds pass; calibration and `FROZEN.json` committed; walks complete
and cross-checked; transported-instance check passes on ≥ 3 curves per source;
the Drive read-back matches. The executor stops and reports here.

### Phase B: measurement (gate: G-B)
1. Run the manifest in order, one `ecbench run --wait` per spec (sessions
   serialize on ecbench's lock; `--cpus auto`; never build or test during
   a session). Order: V3 → V1 → V2 → V4, family D → U → S → N inside each
   volcano, so an early stop still leaves complete strata.
2. After **each** session: `ecbench verify --replay 12` (must be `OK`). Then
   append to `runs/progress.jsonl`. Every 4 sessions, or every 2 h: commit +
   push, upload the new session tarballs to Drive, and read them back.
3. Guards: stop the campaign at 16 h total wall, < 1.5 GiB free disk, or 3
   failed audits in a row. A stopped campaign is a result with its denominator.
4. The transfer-cost microbenchmark (§3) runs once at the end, in its own lock
   window.

**G-B:** every finished session is audited `OK`, pushed and on Drive, and
progress.jsonl's counts match the session directories.

### Phase C: analysis and report (gate: G-C)
1. `ecbench compare` for each (V, src family): `rho-neg` vs each arm, plus
   `rho-neg` vs `rho-neg-aa`. `ecbench table`. Use the §5a pricing band on
   every arm.
2. The statistics in §5, run by `analysis/analyze.py`, written before Phase B
   data are read, from the calibration data only.
3. The report in the research-visuals format: `REPORT.md`, `report.pdf`, the
   volcano diagram coloured by mean `S` per method, cost-by-level forest plots,
   a transfer-cost-vs-solve-cost plot, and a CV/precision table. Copy them to
   `benchmarks/20261017-volcano-ecdlp-census/` and add a row to
   `benchmarks/README.md`.
4. Final upload (git + Drive), then update the PR.

## 5. Statistical design and confidence

**Unit structure.** solve ⊂ target ⊂ (curve, source) ⊂ curve ⊂ level ⊂
volcano. A level's sample is its **curves**, not its solves. Crater levels
with 1–2 curves cannot generalize across curves, so their per-curve variance is
borrowed from the deeper levels, and the report says so.

**Known dispersion (sets n before data).**
- Rho: the collision time is near Rayleigh, CV = √(4/π − 1) ≈ 0.52 (#2030
  measured sd/mean 0.52).
- BSGS textbook: ops = m + j with j uniform on [0, m), so CV ≈ 0.19. With
  shared `k` it is identical across isogenous curves (a sanity check).
- IC: unknown. Estimated in Phase A calibration, and n is re-derived from it
  then (frozen in `FROZEN.json`).

**Precision.** The relative 95% CI half-width of a mean is h ≈ 1.96·CV/√n_eff,
with n_eff = n / (1 + (m−1)·ICC) (m solves per curve; ICC estimated).

| method | CV | n per level for h = 10% | 5% | 3% |
|:--|--:|--:|--:|--:|
| rho | 0.52 | 104 | 416 | 1,155 |
| BSGS | 0.19 | 14 | 56 | 154 |
| IC (CV = 1.0, placeholder) | 1.0 | 384 | 1,537 | 4,268 |

**Targets.** Rho and BSGS are cheap at these sizes. Each level gets ≥ 2,000
rho solves (rounds raised for small levels), so a level-vs-level ratio is
resolved to about ±3.2%. IC gets whatever the 16 h budget allows after Phase A
fixes its cost. The plan reports the h it achieved; it does not claim more.

**Tests (preregistered).**
1. **Level effect, per (volcano, method):** a mixed model
   `log S ~ level + (1 | curve) + (1 | target)`. Report the per-level ratio to
   the crater with a cluster-bootstrap 95% CI: resample curves within a level,
   then solves, 10,000 replicates, seed recorded.
2. **Equivalence (TOST)** at ±5% for rho and BSGS, and at ±10% for IC (or the
   achieved precision if wider, flagged). Report "equivalent", "different" or
   "inconclusive"; never infer "no effect" from p > 0.05.
3. **Curve effect:** the between-curve variance component against the
   isomorphic-control spread (the negative control). A curve is unusual only if
   it sits outside the control's 99% band after Holm correction across curves.
4. **Transfer:** k recovery rate (must be 100%, binomial 95% CI reported).
   Paired D vs U vs N on shared `k_i` (BSGS: an identity check; rho and IC: a
   ratio with CI). Transport cost as a fraction of solve cost.
5. **Positive control:** V3 crater with automorphism-folded rho must reproduce
   #2030's ≈ 1.9× (j = 0) effect. If it does not, the detector is broken and
   the level tests are void.
6. **Multiplicity:** Holm–Bonferroni across (volcano × method × level) for
   tests 1–3. The families are declared here.
7. **Wall time:** reported only when ecbench admits it (`L2`, outside the A/A
   noise). Otherwise it is labelled descriptive. 4-vCPU shared VM: expect L0/L1.
8. **IC vs rho:** per row, the gap to matched `rho-neg` in `S` with CI across
   the 0 / 0.1 / 1 GAE band. That is the "distance remaining", not a pass/fail.

**What would change the conclusion.** A level ratio whose CI excludes 1, that
survives Holm, exceeds the control spread, and replicates on a second volcano
of the same ℓ. One volcano alone gives "observed", not "established".

## 6. Storage

- **git:** `research/volcano_ecdlp_census_20261017/` holds the census (zstd
  JSONL if > 5 MB), specs, session summaries, audits, the analysis, the report
  and the PDF. Raw ecbench session dirs are committed if the total is ≤ 200 MB,
  otherwise tarballs go to Drive with SHA-256 recorded in git.
- **Drive:** the folder `volcano-ecdlp-census-20261017` holds one tarball per
  phase/chunk. `DRIVE.md` lists file id, name, bytes and SHA-256.
- Nothing is deleted or edited after commit. Interrupted sessions stay.

## 7. Compute budget

4 vCPU, 15 GB RAM and about 6 GB of free disk on a shared cloud VM. Phase A
takes ≤ 3 h, Phase B ≤ 16 h wall (guarded), and Phase C ≤ 2 h. Budgets are
advisory. A guard stop is reported with its denominators, never padded.

## 8. Out of scope / not claimed

- No speedup claim. Any apparent IC "win" goes to the breakthrough-tier
  validator before it is mentioned as more than "observed".
- No ecbench `claim` (those are Koblitz-only public-target claims).
- Binary / Koblitz volcanoes; the full all-ℓ class measured; sizes > 26 bits for IC.

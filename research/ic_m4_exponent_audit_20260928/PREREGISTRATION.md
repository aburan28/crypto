# Pre-registration: the `m = 4` exponent audit (algebraic index-calculus decomposition on binary Koblitz curves)

Registered 2026-09-29, before any audit cell ran. Classification: **stage diagnostic**
(AGENTS.md §8). The audit collects no relations, computes no logarithm and gives no end-to-end
ECDLP cost. Its source is `DECOMPOSITION-SURVEY.md` §5
(`research/ic_candidate_tournament_20260915/campaign_20260916/`). The survey fixes the
threshold, the decision rule and the controls. This document fixes everything else before
the data exists.

Phase 1 (this document) did four things: it built the tooling, reproduced the frozen baseline
(§5.1), put one smoke target per arm through the pipeline (§9) and wrote this registration.
**No measured audit cell has run.** The smoke targets used a different seed (7, not
20260928) and are labelled `smoke`. They are not data.

## 1. Object

This is the chained `m = 4` Semaev system of the survey's §3.1:
`S₃(x₁,x₂,e₁), S₃(e₁,x₃,e₂), S₃(e₂,x₄,x_R)` with `S₃(a,c,d) = (a+c)²d² + acd + (ac)² + b`.
It is Weil-descended to a Boolean system with `4ℓ + 2n` unknowns and `3n` cubic equations
(`koblitz_groebner::build_decomposition_system`, `m = 4`).

- **Curves.** `K_a : y² + xy = x³ + a x² + 1` over `F_{2^n}`, `a ∈ {0, 1}`,
  `n ∈ {9, 11, 13, 15, 17, 19}`. The prime sizes are 11, 13, 17 and 19. The anchors are 9
  and 15, the sizes of the frozen `m = 4` cells.
- **Structural exclusion, fixed now.** The survey's grid names three cells that do not exist
  in the tooling. `KoblitzCurve::new` admits a curve only when `#E = h·r` with `r` prime,
  `r² ∤ #E` and `r > h`. That fails for three cells:

  | cell | `#E` | why it fails |
  |:--|:--|:--|
  | `icv1-f2m11-tm67-f393fc83` | 2116 = 2²·23² | `23² ∣ #E` |
  | `icv1-f2m17-t101-e6c4b64d` | 130972 = 2²·137·239 | `r = 239 < h = 548` |
  | `icv1-f2m13-tm181-25e736d2` | 8374 = 2·53·79 | `r = 79 < h = 106` |

  They are not run and not censored: they are absent. **Nine cells remain**, and every `n`
  in the grid has at least one curve:
  `icv1-f2m9-t5-81e744be`, `icv1-f2m13-t181-515ee569`, `icv1-f2m15-tm275-2d22ff5d`, `icv1-f2m19-t797-b6cf2467` and `icv1-f2m9-tm5-4a3ea183`, `icv1-f2m11-t67-05f5aa36`, `icv1-f2m15-t275-b7f03703`, `icv1-f2m17-tm101-00378d4e`, `icv1-f2m19-tm797-9c54981b`.
- **Factor base.** `F = {P : x(P) ∈ V}`, with `V` a uniformly random `ℓ`-dimensional
  `F_2`-subspace of `F_{2^n}` and `ℓ = round(n/4) = ⌊(n+2)/4⌋`, giving 2, 3, 3, 4, 4, 5.
  `V` is drawn exactly as `examples/dreg_ladder.rs` draws it (`random_subspace_basis`,
  reproduced verbatim in the audit binary) from `StdRng::seed_from_u64(cell_seed)`, where
  `cell_seed = 20260928 ^ (a≪48) ^ (n≪40) ^ (ℓ≪32)`. The draw is the first one whose `F` is
  non-empty and `m = 4`-admissible: `m` cofactor classes `[r]P` can sum to `O`.

  Every registered cell took its first draw. The draws were made with `--targets 0` before
  registration, which decides no target. The bases `V` are listed in `cell-objects.txt`.
  Field degree `n`, subgroup order `r`, cofactor `h` and unknown count are kept separate, per
  AGENTS.md §8b:

  | cell | proper intermediate subfields over `GF(2)` | `r` | `h` | `ℓ` | unknowns | `|F|` |
  |:--|:--|--:|--:|--:|--:|--:|
  | `icv1-f2m9-t5-81e744be` | `GF(2^3)` | 127 | 4 | 2 | 26 | 3 |
  | `icv1-f2m13-t181-515ee569` | none | 2,003 | 4 | 3 | 38 | 11 |
  | `icv1-f2m15-tm275-2d22ff5d` | `GF(2^3)`, `GF(2^5)` | 751 | 44 | 4 | 46 | 11 |
  | `icv1-f2m19-t797-b6cf2467` | none | 130,873 | 4 | 5 | 58 | 25 |
  | `icv1-f2m9-tm5-4a3ea183` | `GF(2^3)` | 37 | 14 | 2 | 26 | 1 |
  | `icv1-f2m11-t67-05f5aa36` | none | 991 | 2 | 3 | 34 | 11 |
  | `icv1-f2m15-t275-b7f03703` | `GF(2^3)`, `GF(2^5)` | 211 | 154 | 4 | 46 | 13 |
  | `icv1-f2m17-tm101-00378d4e` | none | 65,587 | 2 | 4 | 50 | 21 |
  | `icv1-f2m19-tm797-9c54981b` | none | 262,543 | 2 | 5 | 58 | 31 |

  The anchors `n = 9, 15` are composite. Their subfield structure and large cofactors are
  disclosed confounds (§10). No arm uses a subfield.

  The Semaev system is `x`-only. Its cost depends on `V` and `x_R`, not on `|F|`. The
  `icv1-f2m9-tm5-4a3ea183` base is a single point, the 2-torsion point: every target there is refuted, but
  the system is still well defined and is still decided.
- **Targets.** 16 per cell, uniform and natural (unplanted): `R = [k]G` with
  `k = 1 + (u mod (r−1))`, where `u` is the first `u64` of
  `StdRng::seed_from_u64(cell_seed·0x9E3779B97F4A7C15 + 0x7A260000 + t)` and `t = 0…15`.
- **Side row not run.** The survey's `n = 31` stable-subspace side row (`ℓ = 5`) has
  `4·5 + 2·31 = 82` unknowns. That exceeds the engine's 64-variable limit (`MAX_VARS`), so it
  cannot be built with this engine. For the same reason `n = 19` is the largest size this
  engine can host at `ℓ ≈ n/4`: `n = 23` would need 70 unknowns.

## 2. Engine

The engine is the inherited F4 with the interleaved chain order and linear elimination of
`RESEARCH_CHAIN_SPLIT_ORDER.md`, which the survey calls "the best in the tree", **as frozen**:

- **Source.** Commit `2809b498f3e45bc0379c0b832c4672c63cae5a73`. That is the commit the frozen
  `research/chain_split_order_20260924` totals were built from (its `manifest.json`), and it
  is the only tree on which they reproduce (§5.1).
- **Added file.** One new file is added to that tree: `examples/m4_exponent_audit.rs` from
  this branch (sha256 `58acf3bd034c101694a81082d4e9acf391a490adb54df6dc2f8e8a69a8f80e20`).
  **No library file is changed**, neither on this branch nor in the frozen tree. The audit
  binary builds `F` by editing the public fields of a library-built base, and copies the
  subspace sampler into the example because it is missing at `2809b498`. The same file also
  builds on the branch head.
- **Lock file.** `Cargo.lock.pinned` (sha256 `469209e8…a869`). `Cargo.lock` is git-ignored in
  the repository.
- **Build.** `WORK=<scratch> research/ic_m4_exponent_audit_20260928/build.sh`, which runs
  `git archive 2809b498` (Cargo.toml, src, examples, benches, docs/ic/calibration.json), adds
  the example, and runs `cargo build --release --locked --example m4_exponent_audit --example
  groebner_stage_bench`. The rebuilt binaries are byte-identical to the phase-1 build:
  - `m4_exponent_audit-2809b498` sha256
    `de2ec1b983b970dc3924cfe459f5bfce726e148dc8cf27d0d4bcdda5480900f2`;
  - `groebner_stage_bench-2809b498` sha256
    `955654969a38086075a3205467a2c05be6eb1d63d56fde2c61b0f2a44bac7fc5`.
- **Engine flags**, all set explicitly on every cell, run under `env -i`:
  - `KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete`;
  - `SolverEngine::default()`, which resolves to `InheritedF4 { max_degree: 3 }` with split
    rule `HighestFree`;
  - reducer `m4ri`;
  - the default size caps (`F4_F2_MAX_ROWS` 20,000 and `F4_F2_MAX_COLS` 40,000) on every
    arm except `degree`;
  - node budget 20,000, as in `groebner_stage_bench`.
- **Oracle.** `koblitz_index_calculus::groebner_decompose(…, m = 4, engine, 20_000)`, the same
  call `groebner_stage_bench` measures. One thread per process.
- **Branch head at registration:** `3ec1b3c7fe24660b75713c79b920655532104b0c` (this
  lane's working branch). While phase 1 was running, a concurrent lane advanced the branch
  from `1722bad1`, where the drift runs in §5.1 were built. `git diff 1722bad1 3ec1b3c7 --
  src examples Cargo.toml` is empty. The audit files are new and not yet committed.

## 3. Arms

| arm | what runs per target | metric | targets per cell |
|:--|:--|:--|--:|
| `semaev` (primary) | `groebner_decompose` on the chained system. Ground truth comes from exhaustive enumeration, uncharged; a returned decomposition is re-added in the group. | 64-bit word XORs in the Macaulay eliminations (`f4_profile().word_ops`, reset per target) | 16 |
| `enumerate` (null 1) | the library's exhaustive `decompose` (non-decreasing index tuples, depth first, stop at the first) with every point addition counted | point additions | 16 |
| `null` (null 2) | `koblitz_bench::random_control_system(n_vars, n_eqs, degree, mean terms/eq, seed_t)`, taking the target's own Semaev system's unknown count, equation count, degree (3) and density. It is solved by the same engine, options and node budget in natural variable order, with every root rejected (a full-tree search, the refutation analogue). `seed_t = cell_seed ^ 0x0C017201·(t+1)` | word XORs | 8 (the first 8 targets) |
| `degree` (secondary) | exact Boolean root count of the chained system (`O(2^{2ℓ+n})`), cross-checked against the engine's own full-tree root count. On the first 4 targets with no root: `solving_degree` in natural layout up to `d_max = 6`, with the `dreg_ladder` caps (`F4_F2_MAX_ROWS=F4_F2_MAX_COLS=50000000`) | refutation degree `D`: resolved, `≥ d_max+1`, or caps | 16 counted, ≤ 4 measured; cells `icv1-f2m9-t5-81e744be`, `icv1-f2m9-tm5-4a3ea183`, `icv1-f2m11-t67-05f5aa36`, `icv1-f2m13-t181-515ee569` (`ℓ = 2, 2, 3, 3`) |

The null has 8 targets rather than 16 for a reason. The smoke target (§9) was censored at
the full 20,000-node budget already at `n = 9`, so the null is expected to censor. Eight
targets decide whether more than half censor at half the CPU. This choice was made
before any audit cell ran.

## 4. Metric

- **Per target:** `T` = word XORs, `semaev` arm.
- **Per cell:** the **lower median**, the ⌈N/2⌉-th smallest of the cell's `T`. A censored
  target counts as `+∞`: it is a lower bound, not a value. The lower median is finite if and
  only if at most half the targets censor, which matches the drop rule in §7.
- **`ĉ`:** the ordinary least-squares slope of `log₂(lower median)` against `n` over every
  retained cell, pooling both curves.
- **Band:** bootstrap over targets, `B = 10,000`, `random.Random(20260928)`. In each
  replicate, every cell's targets are resampled with replacement within the cell, the lower
  medians recomputed, and the slope refitted.
  - A cell whose resampled median is `+∞` is left out of that replicate.
  - A replicate with fewer than 3 distinct `n` is invalid. Invalid replicates are counted.
  - The band is the 2.5th to 97.5th percentile of the valid replicates.
- **Also reported, not decisive:**
  - the refuted-only fit (oracle verdict `refuted`);
  - the satisfiable-only medians. They are fitted only if every retained cell has at least 3
    satisfiable targets. The expected yield is `C(|F|+3,4)/#E`, well under 1 in 16 on most
    cells, so no fit is expected;
  - a prime-`n`-only fit (`icv1-f2m13-t181-515ee569`, `icv1-f2m19-t797-b6cf2467`, `icv1-f2m11-t67-05f5aa36`, `icv1-f2m17-tm101-00378d4e`, `icv1-f2m19-tm797-9c54981b`), read with the same rule. It
    cannot change the verdict. If its reading differs from the primary one, the difference
    is reported.
- **Units.** Word operations are not converted to group operations. The conversion is a
  constant and does not move a slope (survey §5.3).
- **Implementation.** `analyze.py`: pure Python, deterministic, fixed now.

## 5. Controls

### 5.1 Baseline reproduction: PASSED (exact), phase 1

The frozen cells the survey names were re-run through `groebner_stage_bench` built at
`2809b498`, with the `reference` and `candidate` policies, both whole ladders
(`baseline.sh`). Output: `baseline/REPRODUCTION-2809b498.txt` and
`baseline/rerun-2809b498/*/stage.json`.

| cell | arm | `tables.md` (frozen) | rerun | verdict digest |
|:--|:--|--:|--:|:--|
| `icv1-f2m9-t5-81e744be` m=4 (chain) | reference | 2,868,312 | **2,868,312** | `7b29f94a…` = |
| `icv1-f2m9-t5-81e744be` m=4 (chain) | candidate | 917,450 | **917,450** | `85f028df…` = |
| `icv1-f2m15-t275-b7f03703` m=4 (chain-holdout) | reference | 345,384,853 | **345,384,853** | `f859ba30…` = |
| `icv1-f2m15-t275-b7f03703` m=4 (chain-holdout) | candidate | 39,537,587 | **39,537,587** | `edbe5196…` = |

All 14 rows of the two ladders, not only these four, match the frozen rep-1 `stage.json`
exactly in:

- word XORs;
- F4 calls;
- verdict digest;
- decomposed count;
- matrix rows and columns;
- reductions, splits and infeasible branches.

**Disclosed drift: the frozen totals do not reproduce at the branch head.** The same bench
built at `1722bad1`, with the same three policy variables, gives different totals
(`baseline/DRIFT-head-1722bad1.txt`):

| cell | arm | frozen | branch head |
|:--|:--|--:|--:|
| `icv1-f2m9-t5-81e744be` m=4 | reference | 2,868,312 | 595,511 |
| `icv1-f2m9-t5-81e744be` m=4 | candidate | 917,450 | 294,112 |
| `icv1-f2m15-t275-b7f03703` m=4 | reference | 345,384,853 | 83,790,974 |
| `icv1-f2m15-t275-b7f03703` m=4 | candidate | 39,537,587 | 27,620,511 |

- The decomposed counts are the same.
- On `icv1-f2m9-t5-81e744be` the candidate's F4 call count is the same but the verdict digest differs: a
  different decomposition is found.
- The cause: after `2809b498`, 20 commits, merges included, touched the engine sources
  (`koblitz_groebner.rs`, `inherited_f4.rs`, `pq_groebner_f2.rs`, `polynomial_reuse.rs`).
  They include the Four-Russians kernel (`784c0295`), support-local multipliers
  (`1f0d9751`) and set-bit row readback (`0cce6204`), several of them with new defaults.

That is why the registered engine is the frozen tree: it is the only engine on which the
survey's reproduction control can pass. The head engine is 3.1× cheaper than the frozen one
at `n = 9` and 1.43× cheaper at `n = 15` (candidate policy). Its advantage shrinks with `n`
on these two cells, which points to a head slope no lower than the frozen engine's. This is a
two-cell observation, not a fit.

### 5.2 Random-system null (survey §5.5)

The null must grow **faster** than the Semaev systems. Read it at every retained Semaev cell
(one where more than half the Semaev targets are uncensored):

- **pass** if the null's lower median is censored (dearer than its budget) or above the Semaev
  lower median at every such cell;
- **fail** if it is at or below the Semaev lower median at any such cell.

If both arms are measured (uncensored) at 3 or more distinct `n`, their slopes are also
compared, with the `null` bootstrap band reported beside the primary one. Otherwise the
readout says "growth not measurable within budget". That wording is expected, and it is not
a fail.

**Gating:**

- A **fail** downgrades an *alive* verdict to *inconclusive (null control failed)*.
- A *closed* verdict does not rest on structure and stands.

### 5.3 Enumeration null (survey §5.5)

The survey's `c = (m−1)x ≈ 0.75` is the asymptotic value. For refuted targets on the
registered bases, full-tree point additions are a function of `|F|` alone, and their LS slope
on `n` is **0.828** (`|F|` from §1).

- **pass** if `ĉ_enum ∈ [0.60, 1.00]`.
- Outside that range, the target and subspace generation is suspect. Every verdict is then
  reported as *instrument check failed* and nothing is promoted.

### 5.4 Instrument checks (fixed now)

- **Oracle verdict against enumeration truth.** Every disagreement is listed. The split uses
  the oracle's verdict, since it is the oracle's cost path being measured. A disagreement
  shows an oracle that misses degenerate decompositions. It is not a measurement error.
- **Group verification.** A returned decomposition that fails the group re-addition **voids
  the audit**.
- **Degree arm count check.** A mismatch between the exact count and the engine's complete
  full-tree root count voids the `degree` arm only.

## 6. Decision rule (verbatim from survey §5.4)

- **alive** if the band's upper end is below `c* = 0.25`;
- **closed for this engine at these sizes** if its lower end is above `0.25`;
- otherwise **inconclusive**.

This is applied to the primary fit (§4), after the gating in §5.2 to §5.4. If fewer than 3
distinct `n` are retained, the fit is not made and the verdict is **inconclusive
(insufficient uncensored sizes)**. `analyze.py` applies the gating and prints the result as
`REGISTERED VERDICT`.

**Secondary check (non-decisive):** the refutation degree `D(ℓ)` on the degree cells, with
`s* ≈ 0.05` (survey §2).

- `s` is the LS slope of resolved refutation degrees on `ℓ`.
- Lower bounds (`≥ d_max+1`) and cap hits are listed, not fitted.
- With only `ℓ ∈ {2, 3}`, any `s` is a two-level read.
- The smoke target (§9) already reached `D ≥ 6` at `ℓ = 2`, so lower bounds are the likely
  outcome.

## 7. Stop rules and censoring

- **Node budget (deterministic).** 20,000 per target for `semaev` and `null`. Hitting it gives
  `exhausted` and the target is **censored**. The engine's size caps (`oversize`) are part of
  the engine and are recorded per target.
- **CPU-seconds limit per cell (`ulimit -t`; machine protection, not a scientific rule):**

  | arm | limit per cell |
  |:--|:--|
  | `semaev` | 120 s (`n ≤ 15`); 300 s (`n = 17`); 600 s (`n = 19`) |
  | `enumerate` | 60 s |
  | `null` | 120 s |
  | `degree` | 300 s, plus `ulimit -v` of 10 GB |

  A killed process keeps every line it wrote. The targets it never reached are **censored**.
- **Censoring is never negative evidence.** A censored target is a lower bound. A cell that
  censors on **more than half** its targets **drops out of the fit**, and the drop is
  recorded in the readout. The same rule applies to both nulls.
- No cell is re-run with a larger budget, and no target is replaced, under this
  registration. Any change is a dated, additive amendment below this document. Nothing here
  is rewritten.
- Outputs are write-once: the binary refuses an existing `--out`, and `run_audit.sh` refuses
  an existing runs directory.

## 8. Phase 2: exact commands, seeds and CPU budget

**Seeds.**

| seed | value |
|:--|:--|
| audit | 20260928 (all cells, all arms) |
| bootstrap | 20260928 |
| smoke (phase 1 only, not data) | 7 |

**Commands**, from the repository root:

```sh
export WORK=<phase-1 scratch directory>   # holds target/ (reused), src2809/ and bin/
research/ic_m4_exponent_audit_20260928/build.sh           # ~2 min; prints four sha256s, which must equal §2's
# (if a binary hash differs, record both in a dated amendment before any cell runs)
export BIN=$WORK/bin/m4_exponent_audit-2809b498
research/ic_m4_exponent_audit_20260928/run_audit.sh research/ic_m4_exponent_audit_20260928/runs
# = the 31 cell commands of phase2-commands.txt (27 four-way parallel, then 4 serial),
#   then: python3 research/ic_m4_exponent_audit_20260928/analyze.py research/ic_m4_exponent_audit_20260928/runs
```

`phase2-commands.txt` is the exact output of `run_audit.sh --print` and lists every cell
command. Their form is:

```sh
( ulimit -t CPU; exec env -i PATH="$PATH" KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete \
    $BIN --arm ARM --a A --n N --targets T --first 0 --seed 20260928 --label audit [--node-budget 20000] \
    --out runs/ARM/K{A}n{N}l{ELL}.jsonl ) 2> runs/ARM/CELL.stderr; echo $? > runs/ARM/CELL.exit
```

The `degree` cells add `ulimit -v 10000000`, `F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000`,
`--d-max 6 --max-unsat 4`.

**CPU (advisory).** There are 27 cells in the parallel set: 9 `semaev`, 9 `enumerate` and
9 `null`. There are 4 `degree` cells.

- **Hard ceiling (the sum of the per-cell limits):**

  | arm | ceiling |
  |:--|--:|
  | `semaev` | 2,220 s |
  | `enumerate` | 540 s |
  | `null` | 1,080 s |
  | `degree` | 1,200 s |
  | **total** | **5,040 s ≈ 84 CPU-minutes** |

- **Expected, a model and not a measurement:**
  - `semaev`: about 1–2 min, if the frozen readout's ≈0.8 bits per unit `n` holds. That puts
    `n = 19` at about 1 s per target, from 0.115 s per target on the frozen `icv1-f2m15-t275-b7f03703` m=4
    cell. At 2 bits per unit `n` the `n = 19` cells run into their limits, about 20 min.
  - `null`: about 10–18 min. Most cells are expected to hit the 120 s limit or censor at
    20,000 nodes; the smoke target used 1.2 s at `n = 9`.
  - `degree`: about 5–20 min.
  - `enumerate`: under 1 s.
  - **Expected total: about 20–45 CPU-minutes.** Wall time is about 10–25 min with 4 lanes
    plus the serial degree cells.
- **Disk:** under 5 MB of JSON lines. The build reuses `$WORK/target` (≈250 MB).

## 9. Smoke (phase 1; labelled `smoke`; not data)

One target through each arm, at `icv1-f2m9-t5-81e744be` with seed 7: `ℓ = 2`, `V = [113, 44]`, `|F| = 3`,
target `k = 24`, 26 unknowns, 27 equations, degree 3, 54 terms per equation. Files are in
`smoke/`.

| arm | outcome |
|:--|:--|
| `semaev` | refuted; enumeration truth: not decomposable (agree). 72,244 word XORs, 20 F4 calls, 5 splits, 0.005 s |
| `enumerate` | refuted; 19 point additions |
| `null` | **censored** at the 20,000-node budget: ≥ 15,772,537 word XORs, 19,999 F4 calls, 1 root seen (the planted all-zero root), 1.17 s |
| `degree` (default caps) | 0 Boolean roots (exact count = engine full-tree count = 0). The caps stopped it at degree 4 (`caps_hit`) |
| `degree` (dreg caps, `d_max = 5`, same target re-run once) | `at_least 6`, 1.0 s |

The second `degree` line is why the `degree` arm registers the `dreg_ladder` caps and
`d_max = 6`.

## 10. What it cannot show

(Survey §5.7, plus what this registration adds.)

- **No asymptotic statement.** The claim at stake is "∃ solver family, ∀ large `n`,
  `T ≤ 2^{cn}` with `c < 1/4`", and finite data can only falsify within the tested family and
  sizes (survey §3.1, quantifier order).
- **Nothing about another solver family,** or about the branch-head engine, whose totals
  differ (§5.1).
- **Nothing at `n ≥ 131`.** Any later ECC2K-130 transfer claim owes the `m = 83` gate of
  AGENTS.md §8a.
- **An *alive* verdict is a reason to scale, not a result.**
- **Nothing beyond `n = 19` or `m = 4` for this engine,** which is capped at 64 unknowns.
- **The fit mixes several things:**
  - two curves;
  - composite and prime `n`: the anchors 9 and 15 carry subfield structure and cofactors of
    4–154;
  - a jittered ratio `x = ℓ/n ∈ [0.22, 0.27]`;
  - bases with 1 to 31 points.

  None of these moves the x-only system's cost directly, but each is a confound on a 9-point
  fit.
- **A "closed" verdict closes `m = 4` for this engine at these sizes. It does not close the
  route** (survey §0.5, §5.6).

## 11. Successor and revisit conditions

These come from survey §5 and are kept as written:

- **`m = 5`** (`c* = 0.30`) is the next rung. The chain adds `n` unknowns per summand, which
  makes `γ*` smaller.
- **Any engine that measures below 0.25 at `m = 4`** reopens this.
- **A published degree bound** for Weil-descended chains, in either direction, reopens it.
- **§3.2's symmetry lever** reopens it, if it lowers the slope `s` rather than cutting `D`
  by a constant.

Added by this registration:

- **Branch-head engine.** A rerun of these nine cells on the branch-head engine is a cheap
  successor if the verdict is inconclusive, or to test whether the post-freeze changes move
  the slope. It needs its own reproduction baseline, which would be the head totals in
  `baseline/head-1722bad1`.
- **More than 64 unknowns.** An engine beyond 64 unknowns is needed for `n ≥ 23` at
  `ℓ ≈ n/4`, and for the survey's `n = 31` stable-subspace side row. It is the precondition
  for turning an *inconclusive* result into a verdict by adding sizes.

## Files

All paths are relative to `research/ic_m4_exponent_audit_20260928/` unless stated.

| file | what it is |
|:--|:--|
| `examples/m4_exponent_audit.rs` (repository root; new) | the audit binary, all four arms |
| `build.sh` | builds the registered engine |
| `baseline.sh` | the reproduction control |
| `run_audit.sh` | phase 2 |
| `analyze.py` | the registered readout |
| `phase2-commands.txt` | the exact phase-2 command list |
| `cell-objects.txt` | the registered `V` and `|F|` for every cell |
| `Cargo.lock.pinned` | the pinned lock file |
| `baseline/` | the reproduction outputs and the head drift |
| `smoke/` | the smoke outputs, not data |

## Amendments

(None. Amendments are appended here, dated, and never edit the sections above.)

### Amendment 1 — 2026-09-29, before any audit cell ran: `rustfmt`, and the binary hash

The repository's `rustfmt (changed files)` CI check failed on
`examples/m4_exponent_audit.rs`. The file was reformatted with `rustfmt --edition 2021`.
Only whitespace and line breaks changed; no token did.

- **Source sha256:** `58acf3bd…0e20` became `c13e9efc8579c648c32ce5d6703131655e64e156dd7197e078dcca8a2e0bfdf8`.
- **`m4_exponent_audit-2809b498` sha256:** `de2ec1b9…00f2` became
  `79ccfdf395bc22776ce061921851386dc47038a7d1b0fe3a7dda44e139acb55e`. The difference is the
  line numbers compiled into panic locations.
- **Unchanged:** `groebner_stage_bench-2809b498` (`95565496…7fc5`) and `Cargo.lock.pinned`
  (`469209e8…a869`).
- **Behaviour checked equal.** The new binary re-ran the §9 smoke targets for `semaev`,
  `enumerate` and `null` (seed 7, `icv1-f2m9-t5-81e744be`). Every field of every report equals the
  committed `smoke/` output except `secs`, the wall time. That includes word operations, F4
  calls, matrix rows and columns, reductions, splits, verdicts and the null's censoring
  point.

The binary with the new hash is the registered engine. Nothing else changes.

### Amendment 2 — 2026-09-29, before any audit cell ran: `clippy` on the module doc comment

The repository's `clippy (crypto)` CI check (`cargo clippy --all-targets -- -D warnings`,
Rust 1.98) failed with 18 `doc_overindented_list_items` errors in the module doc comment of
`examples/m4_exponent_audit.rs`. The continuation lines of the four arm descriptions are
re-indented from 17 spaces to 2. Only comment whitespace changes; no line is added or
removed.

- **Source sha256:** `c13e9efc…bfdf8` (Amendment 1) became
  `5bfe99f0ece7d8fe1345b947136948b7588bd1179a02aa8a7770b744ffab1cc8`.
- **The binary is unchanged.** `m4_exponent_audit-2809b498` rebuilt by `build.sh` is
  byte-identical to Amendment 1's
  (`79ccfdf395bc22776ce061921851386dc47038a7d1b0fe3a7dda44e139acb55e`), and so is
  `groebner_stage_bench-2809b498` (`95565496…7fc5`). The registered engine is unchanged.
- **Checked locally** with clippy 0.1.94 at `-D warnings`: the previous file fails with the
  same lint, and this one passes.

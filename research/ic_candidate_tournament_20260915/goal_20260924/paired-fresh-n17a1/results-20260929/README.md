# Retained n17a1 paired correctness panel

Classification: **accounting and correctness diagnostic**. All six local arms
are terminal, including the separate SAT v1 control. No arm is eligible for
promotion and every admissible comparative online cost, speedup, normalized
`S`, operation ratio to rho and operation ratio to a floor remains **unknown**.
The three earlier sealed confirmation rounds remain closed.

The allocated public target is `[52411,72106]` on
`EC1N17Ckb1hbbe2b5b6b1e6`, over `GF(2^17)` with subgroup order `65587`.
The field degree is not a subgroup bit count. Fixture seed `2026092948` is
provenance only; the workers received a point and no secret scalar. SAT v2,
the local pair-table incumbent and both rho diagnostics independently recover
`24886` and replay `[24886]G = Q`. The F5 arm finishes its registered collection
budget without full rank and never enters target descent.

## One-target result table

The raw online observations below use **nanoseconds throughout**. They are
retained instrument readings on a contended macOS ARM64 host, not calibrated
comparative costs. In completed arms the online interval excludes reusable
preparation and ends immediately after scalar replay. The independent external
Python audits are separate. Unknown values are not zeros.

| Arm | Actual usable base / folded columns | Ordinary queries / verified relations / rank | Target attempts | Verified scalar | Raw online observation, ns | Admissible online cost, ns | Rho / IC speedup |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Static SAT v2 | 62 / 29 | 189 / 57 / 29 | 2 | 24886 | 72,019,683,375 | unknown | unknown |
| Bounded MatrixF5 | 62 / 29 | 256 / 61 / 28 | 0 | none | unknown | unknown | unknown |
| Accepted-source pair-table IC | 272 / 8 | 9 / 9 / 8 | 1 | 24886 | 25,791 | unknown | unknown |
| Selected-source signed-Frobenius rho | not applicable | not applicable | one target | 24886 | 158,125 | unknown | not applicable |
| Generic rho diagnostic | not applicable | not applicable | one target | 24886 | 142,750 | unknown | not applicable |

The [machine-readable table](DIAGNOSTIC-PANEL.json) retains all row keys,
exclusive online phases, process observations and null claim fields. Its raw
observations come from the untouched audits and result files in the archive.
It is a post-execution analysis, not a new candidate or measurement registration.

SAT v2 and F5 share the dimension-six standard-subspace geometric construction:
63 geometric points map to 62 distinct nonidentity subgroup points, with
29 sign/Frobenius columns. Their usable-set digest is
`3461ea31b5bb93065826e34b74658ba916477971c42a212e114387075eb4a2d6`.
The incumbent changes both algorithm and factor-base policy. Its nominal bound
102 produced 272 actual usable points and eight columns; sampler overshoot is
retained in the identity rather than mislabeled as `fb102`.

| Arm | Candidate or reference ID | Workload ID | Run number |
| --- | --- | --- | --- |
| SAT v2 | `IC1N17Ckb1fb62PDP3satRCsampleLAgaussTDpdpISO0h55d1911e6064` | `70fee57abbcf` | 0 |
| F5 | `IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h57ee76cbeaed` | `4937b2d124b2` | 0 |
| Incumbent | `IC1N17Ckb1fb272PDP3pairRCwalkLAgaussTDpdpISO0hb6019b409f77` | `ee5afbda6851` | 0 |
| Selected rho | `RHO1N17Ckb1h90e231acc61f` | `e0350dc3e3a4` | 0 |
| Generic rho | `RHO1N17Ckb1hf921d8edfd97` | `98f1c920f1dd` | 0 |

Each complete run ID is the candidate/reference label plus `W<workload-id>R0`;
the exact immutable records are retained. Same point and curve do not silently
make differing cache/resource declarations the same workload identity.

## Natural attempts and the missing F5 direction

SAT v2 retains 57 verified witnesses, 125 source-reported UNSAT attempts and
seven conflict-budget-inconclusive attempts. The independent full-point oracle
checks feasibility for each input; exactly 57 collection inputs are feasible.
It verifies all witnesses, final matrix, base logs, two target attempts and
scalar replay. The source UNSAT labels are not separate checked SAT proof
certificates. All 57 matrix rows are retained; 28 are dependent and none are
duplicates. The observed witness fraction's descriptive Wilson 95% interval
is `[0.240644,0.370436]`, from the original audit. Adaptive rank stopping and
one query stream prevent interpreting this as a universal natural-yield rate.

F5 retains 61 witnessed attempts and 195 bounded-incomplete attempts across
all 256 ordinary queries, with no proved-negative classifications. Its
registered degree-three Macaulay engine uses an F5 row criterion, node budget
4096, batches of eight and one Rayon worker; it is not a complete incremental
Gröbner-basis algorithm. Independent reconstruction verifies every admitted
relation, query, batch and matrix. Sixteen final relation-LA attempts reach
rank 28/29, last increasing at query 236. The zero-descent, exit-code-two result
is `AUDITED_BOUNDED_INCOMPLETE`; its missing recovery phases prevent a
complete-solve cost.

The first 189 F5 queries exactly match SAT's query scalars and relation columns:

| F5 outcome | SAT audited outcome | Matched inputs |
| --- | --- | --- |
| witnessed | valid point witness | 45 |
| bounded incomplete | valid point witness | 12 |
| bounded incomplete | source UNSAT, oracle infeasible | 125 |
| bounded incomplete | conflict-budget inconclusive | 7 |

This prefix ends at SAT's full-rank stop. It is a **post-hoc matched-input
diagnostic**, not a preregistered rate comparison or a selection gate.
The twelve SAT-valid/F5-incomplete query numbers are
`7,13,45,48,50,86,92,99,113,164,173,188`.

A separate modular elimination reconstructs rank 28 over the prime scalar
field. Its nullspace is the unit vector at zero-based column 27, whose point
representative is `[18941,32935]`. Every retained F5 row has zero coefficient
in that column. SAT's rows from queries 164 and 173 have nullspace dot products
51777 and 13810 respectively and cover that missing direction. Both queries
were budget-incomplete in F5. This locates the observed shortfall in PDP
coverage before final LA; it does not support calling the matrix kernel slow
or incorrect. No measured query is retried to replace the failed registration.

## Timing, source and reference limits

- SAT executed before F5, contrary to the allocated order. Its workload declares
  warm preparation while the other records declare cold preparation. Preserve
  these differences; do not collapse their workload IDs or reconstruct a
  favorable paired order after seeing outcomes.
- Host-load sidecars retain loads around 40–54 on 14 logical CPUs. Source,
  test and Git activity also overlapped the long runs. The repository's Linux
  CPU-affinity, `/proc` and PSI isolation gates were unavailable. There is no
  local A/A calibration or repeated-pair uncertainty estimate. In particular,
  the short native intervals cannot substantiate a speedup ratio.
- Original SAT v1/v2 Python manifests omitted the imported namespace-package
  files `producer/evidence.py` and `producer/timing.py`. Their bytes are recovered
  from the recorded measured commits, **after execution**, in the separate
  [historical runtime bundle](../../static-sat-runtime/README.md). Exact replay
  does not retroactively establish complete preexecution source coverage.
  No original manifest, seal, candidate ID or result is rewritten.
- The incumbent and selected rho retain the accepted pair-table source manifest,
  complete local dependency inventories and macOS release binary SHA-256
  `b3d8afacc746a34cf15562452dbec3d75deadbb9faf16051e97ad8fa07cbbd66`.
  This is a separately built portable local binary, not the qualified Linux
  instruction-count binary. The generic worker's source and build binding
  remain separate and replay through its existing admission checker.
- Selected rho requests four interleaved walks, independently checked to clip
  to one effective walk at this group size. It records 19 iterations and one
  restart. Interleaved width is not a CPU-worker count; all local workers use
  one Rayon thread. Collision-table allocation memory is unknown. Retained
  whole-process high-water RSS does not fill that gap.
- F5's original external audit took 78,188,583 ns. The other original external
  Python audit durations were not recorded and remain null. Whole-panel replay
  has its own separately recorded duration. No unrecorded audit cost is inserted
  into an original online interval or treated as zero.
- This one toy point establishes no population-level speedup, common-operation
  `S`, exponent fit, global optimum, m=83 gate or ECC2K-130 improvement.

The separate SAT v1 control solved `[114119,85674]` as scalar 61885 after
131 ordinary attempts: 31 witnesses, 93 source UNSAT and seven independently
reclassified conflict-budget-inconclusive outputs. Its original endpoint also
included a final progress write after scalar replay, so its 3,045,557,458 ns
observation is overinclusive and correctness-only. It is not a paired row for
the fresh point. All its raw files and original audit are retained.

## Archive and independent replay

The committed [evidence.tar.gz](evidence.tar.gz) is 28,948,695 bytes with
SHA-256 `039e535fa73e5c794aeff666bc02989076f04abbf751b5f609884a2785760678`.
[receipt.json](receipt.json) inventories all 3,851 retained regular files.
It includes every closed arm file, failed/inconclusive attempt, native binary,
registration, build, preparation receipt and host-load sidecar, plus the
accepted incumbent source bytes. The Git version of this directory is the
durable location; scratch paths inside historical receipts are provenance.

The three [independent SAT auditor sources](audit-source) are separately
recovered from commit `6958067b53f110be39a7cbcd9627aab43908f906`, with hashes
and recovery scope in [auditor-receipt.json](auditor-receipt.json). Recovery
does not change original candidate identity. The panel replayer uses the
frozen historical SAT runtime, retained raw output and these auditor bytes;
it hashes the Mac workers but never executes a measured producer or SAT solver.

From the repository root:

```sh
python3.12 research/ic_candidate_tournament_20260915/replay_paired_n17_evidence.py \
  --out /private/tmp/ic-paired-n17-new-replay.json
PYTHONPATH=research/ic_candidate_tournament_20260915 python3.12 -m unittest \
  test_paired_n17_evidence
tar -xzf research/ic_candidate_tournament_20260915/goal_20260924/paired-fresh-n17a1/results-20260929/evidence.tar.gz \
  -C /absolute/path/to/new-empty-directory
```

Use a new replay output path. [REPLAY-v2.json](REPLAY-v2.json) is the stronger
final receipt; the earlier [REPLAY.json](REPLAY.json) remains preserved.
Original SAT audit fields still say `same_point_rho_audited: false` because
rho had not then completed. The new joint receipt says true after checking
the retained same-point references. That resolves the reference-audit gap,
while source coverage, host calibration and ordering still block admission.

## Decision and continuation

Do not promote this panel, overwrite a run or redispatch its registrations.
Keep `[52411,72106]` and the SAT v1 point excluded from fresh qualification.
The [next SAT source protocol](../../static-sat-runtime-v3/PROTOCOL.md) requires
a package-complete preexecution archive and imported-module gates in a
versioned registrar, runner and auditor. F5 needs a separately frozen coverage
or budget policy that can reach full rank on disclosed development controls.
Then register new excluded public targets, exact arm limits, the strong
same-point reference and calibrated Linux execution with A/A and declared
order. Complete one-target solves and failed-attempt accounting for both
families are still pending; the persistent research goal remains active.

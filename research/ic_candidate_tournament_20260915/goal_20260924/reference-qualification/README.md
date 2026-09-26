# Archived reference qualification, September 25

The declared development panel completed: **1,290/1,290 native/profile pairs**,
including 30 A/A controls, 315 smoke pairs and 945 development pairs. All 21
configurations completed 45 development jobs each (15 distinct targets, three
process repetitions). No measured timeout, OOM or failed answer occurred in this
panel. The archived control history separately retains three intentional native
failures. This is reference qualification: **zero improvement rounds used, no
promotion, no confirmation data generated or used**.

Implementation [PR 761](https://github.com/aburan28/crypto/pull/761) merged at
`e9b540eef8775f8d0d65e24c319174d8fefa7460`. The measured head was
`3b1112b65f11c6d3678e8bea26c0a22298e2efa1`; its implementation files match the
merge. The [completed run](https://github.com/aburan28/crypto/actions/runs/36127866931)
executed the previously frozen [protocol](PROTOCOL.md). Measurement used Linux
amd64, CPU 3 affinity, one worker thread, 8 GiB address-space cap, 180-second child
limits, Rust 1.94.1 and Valgrind 3.22.0. OS caches were uncontrolled. This panel
contains five exact curve/subgroup cells, not five monotonically increasing
subgroup sizes. Exact fields, groups, public points, canonical identities,
source/build hashes and exclusive stages are retained in the archive.

## Selection and limits

The frozen selection rule chooses `pairinv` for both IC cold instructions and
online native time. It chooses `rho_incumbent_4` for rho cold instructions and
`rho_pairinv_4` for rho online time; retain both references. Requested width 4
executes widths **1, 1, 4, 4, 2** on cells n17a1, n19a0, n23a0, n23a1 and n31a0.
Requests 8, 16 and 32 all clip to **1, 1, 4, 5, 2**. They are separate measured
configurations, not independent algorithms or additional target samples.

`pairinv`'s online ratio to the old `both` incumbent is 0.9781, with descriptive
95% interval [0.9249, 1.0302]. Its complete instruction ratio is 0.9845, while its
cold native ratio is 1.0064. These observations do not meet the 20% goal, and the
online interval includes no improvement. The selected online rho reference costs
8.5202 times `pairinv` online on this panel; that ratio excludes reusable IC/rho
preparation, includes target-dependent work through scalar replay, and is not a
cold-start gain or a confirmed comparison. Cold native times remain close to one
millisecond and include launch/input/report overhead. All rankings are specific
to these archived implementations, development points and this host.

The current general `icx` entry point is not admitted to this public-target
panel. F4/F5/SAT adapters and new factor-base families were not measured here.
The row kernel was checked against independent scalar-field RREF/rank/solution
controls; its full-row arithmetic was not optimized in this qualification.
The allocator over-alignment correction is correctness work. No result here
establishes a globally optimal implementation or a mathematical advance.

## Independent replay and transport finding

All 1,290 receipts pass independent certificate/accounting replay on macOS with
Python 3.12. Exact floating-summary replay there rejects two values: one smoke
bootstrap endpoint and one per-cell geometric mean differ by one ULP each
(relative errors below 1.3e-16). No integer, source, input, completion, identity,
selection or correctness value differs. The diagnostic and its original script
are retained in `audits/` inside the archive. This is a cross-platform numerical
reproducibility limitation, not a pass of the strict whole-bundle audit.

The dedicated Linux evidence-replay workflow performs fresh extraction and all
**22 strict frozen audits**, including exact summary/selection comparisons,
then reproduces every exported row. Its result must pass before this evidence PR
can merge. The original bytes and exact comparisons remain unchanged.

## Primary single-target online measurements

Numbers below are geometric means of per-point medians, with equal weight per
cell. Each job solves one public point; these are not batch-amortized timings.
The aliases denote complete source/configuration policies across cells, not
canonical candidate IDs. [RESULTS.json](RESULTS.json) maps them to the per-cell
IC1/RHO1 identities and retains all 945 development run records with canonical
workload/run IDs. Every row has 45/45 verified development jobs. Intervals are
descriptive development intervals; selection does not provide a familywise
confidence claim.

| Configuration alias | Online ms | / old IC | Descriptive 95% CI | Selected rho / online | Verified |
|---|---|---|---|---|---|
| incumbent | 0.03527943 | 1 | ['1', '1'] | 8.333291 | 45 |
| scaled | 0.03658234 | 1.036931 | ['0.9816509', '1.107543'] | 8.036496 | 45 |
| pairinv | 0.0345057 | 0.9780685 | ['0.924907', '1.030152'] | 8.520151 | 45 |
| rho_incumbent_1 | 0.344969 | 9.778191 | ['7.725357', '12.07888'] | 0.8522324 | 45 |
| rho_incumbent_2 | 0.3160875 | 8.959541 | ['7.392714', '10.75512'] | 0.9301025 | 45 |
| rho_incumbent_4 | 0.2966635 | 8.408965 | ['7.10053', '9.967899'] | 0.9910008 | 45 |
| rho_incumbent_8 | 0.3011474 | 8.536061 | ['7.073196', '10.23016'] | 0.9762455 | 45 |
| rho_incumbent_16 | 0.3032368 | 8.595285 | ['7.11568', '10.25388'] | 0.9695189 | 45 |
| rho_incumbent_32 | 0.3034023 | 8.599977 | ['7.231951', '10.20722'] | 0.9689899 | 45 |
| rho_scaled_1 | 0.3460468 | 9.808741 | ['7.717053', '12.19424'] | 0.8495781 | 45 |
| rho_scaled_2 | 0.3201308 | 9.074149 | ['7.427558', '11.00426'] | 0.9183551 | 45 |
| rho_scaled_4 | 0.3022875 | 8.568378 | ['7.282983', '9.996674'] | 0.9725635 | 45 |
| rho_scaled_8 | 0.2988046 | 8.469655 | ['7.118939', '10.09077'] | 0.9838997 | 45 |
| rho_scaled_16 | 0.306974 | 8.701216 | ['7.395277', '10.30697'] | 0.9577157 | 45 |
| rho_scaled_32 | 0.3058821 | 8.670267 | ['7.259575', '10.32453'] | 0.9611344 | 45 |
| rho_pairinv_1 | 0.3419798 | 9.693462 | ['7.649502', '12.07019'] | 0.8596817 | 45 |
| rho_pairinv_2 | 0.3222622 | 9.134563 | ['7.601432', '10.8611'] | 0.9122814 | 45 |
| rho_pairinv_4 | 0.2939938 | 8.333291 | ['7.066303', '9.854115'] | 1 | 45 |
| rho_pairinv_8 | 0.3024054 | 8.571718 | ['7.197861', '10.14227'] | 0.9721844 | 45 |
| rho_pairinv_16 | 0.3044545 | 8.629802 | ['7.275892', '10.21295'] | 0.965641 | 45 |
| rho_pairinv_32 | 0.2981093 | 8.449944 | ['7.124324', '9.997048'] | 0.9861948 | 45 |

The following paired rows show the chosen IC and online-rho settings on every
individual public target. Full public coordinates, candidate/reference IDs,
workload IDs, all six run IDs and the timing boundary are in
`RESULTS.json:single_target_online`; all three IC sources are retained there.

| Public-target case | IC online ms | rho online ms | rho / IC | Verified |
|---|---|---|---|---|
| n17a1-000 | 0.027832 | 0.243165 | 8.73689 | yes |
| n17a1-001 | 0.025648 | 0.209192 | 8.15627 | yes |
| n17a1-002 | 0.029796 | 0.179415 | 6.02145 | yes |
| n19a0-000 | 0.029495 | 0.230742 | 7.82309 | yes |
| n19a0-001 | 0.030217 | 0.212819 | 7.04302 | yes |
| n19a0-002 | 0.025628 | 0.207037 | 8.07855 | yes |
| n23a0-000 | 0.035396 | 0.491149 | 13.87583 | yes |
| n23a0-001 | 0.032070 | 0.281347 | 8.77290 | yes |
| n23a0-002 | 0.032651 | 0.305421 | 9.35411 | yes |
| n23a1-000 | 0.037911 | 0.259145 | 6.83561 | yes |
| n23a1-001 | 0.036208 | 0.496288 | 13.70658 | yes |
| n23a1-002 | 0.086402 | 0.364382 | 4.21729 | yes |
| n31a0-000 | 0.037860 | 0.411249 | 10.86236 | yes |
| n31a0-001 | 0.042650 | 0.468917 | 10.99454 | yes |
| n31a0-002 | 0.033943 | 0.299360 | 8.81949 | yes |

## Supplementary complete cold native time

Cold means the whole worker process, with reusable preparation and external
launch/input/report/exit overhead charged. Fixture construction remains outside
both algorithms. This boundary is separate from the primary online interval.

| Configuration alias | Cold ms | / old IC | Descriptive 95% CI |
|---|---|---|---|
| incumbent | 0.9548849 | 1 | ['1', '1'] |
| scaled | 0.974083 | 1.020105 | ['0.9906293', '1.052218'] |
| pairinv | 0.9610005 | 1.006405 | ['0.9651924', '1.049098'] |
| rho_incumbent_1 | 1.004023 | 1.05146 | ['1.006284', '1.104561'] |
| rho_incumbent_2 | 0.9697427 | 1.01556 | ['0.9749967', '1.058701'] |
| rho_incumbent_4 | 0.9405923 | 0.9850321 | ['0.9303436', '1.036712'] |
| rho_incumbent_8 | 0.9478406 | 0.9926229 | ['0.9463485', '1.039489'] |
| rho_incumbent_16 | 0.9503984 | 0.9953015 | ['0.9505123', '1.040278'] |
| rho_incumbent_32 | 0.9475625 | 0.9923317 | ['0.9479119', '1.036316'] |
| rho_scaled_1 | 1.01234 | 1.06017 | ['1.024653', '1.105298'] |
| rho_scaled_2 | 0.9826127 | 1.029038 | ['0.9887804', '1.069209'] |
| rho_scaled_4 | 0.9492754 | 0.9941255 | ['0.9456095', '1.042069'] |
| rho_scaled_8 | 0.9506862 | 0.9956029 | ['0.9444905', '1.044679'] |
| rho_scaled_16 | 0.9611215 | 1.006531 | ['0.9565629', '1.058752'] |
| rho_scaled_32 | 0.9497731 | 0.9946467 | ['0.9465444', '1.044426'] |
| rho_pairinv_1 | 1.007289 | 1.05488 | ['1.011843', '1.107358'] |
| rho_pairinv_2 | 0.9787077 | 1.024948 | ['0.9908894', '1.059673'] |
| rho_pairinv_4 | 0.9536487 | 0.9987054 | ['0.9408598', '1.052529'] |
| rho_pairinv_8 | 0.9553195 | 1.000455 | ['0.9534799', '1.046117'] |
| rho_pairinv_16 | 0.9510661 | 0.9960008 | ['0.9515908', '1.039741'] |
| rho_pairinv_32 | 0.9488646 | 0.9936953 | ['0.9536563', '1.035123'] |

## Supplementary complete cold instructions and boundaries

The unit is `valgrind-3.22-amd64-Ir`, not curve additions or field operations.
`S` therefore must not be compared numerically to a rho constant quoted in a
different operation unit. The K-instruction floor is the predeclared weak bound
for this full-rank IC collector; it is inapplicable to rho and cannot establish a
generic lower bound or an algorithmic advance. Instruction costs include all
exclusive phases and unsuccessful target-dependent attempts. Raw stage records
retain base census, query status mix, novel rank, matrix construction/reduction,
descent and scalar replay. Missing costs were not replaced with zero.

| Configuration alias | Cold Ir | S = Ir / sqrt(r) | / old IC | / selected rho | / K floor |
|---|---|---|---|---|---|
| incumbent | 1780785 | 2223.299 | 1 | 0.7138638 | 222598.2 |
| scaled | 1780670 | 2223.155 | 0.9999353 | 0.7138176 | 222583.8 |
| pairinv | 1753224 | 2188.889 | 0.9845231 | 0.7028154 | 219153 |
| rho_incumbent_1 | 2967440 | 3704.83 | 1.666366 | 1.189558 | not applicable |
| rho_incumbent_2 | 2728287 | 3406.248 | 1.532069 | 1.093689 | not applicable |
| rho_incumbent_4 | 2494573 | 3114.458 | 1.400827 | 1 | not applicable |
| rho_incumbent_8 | 2540371 | 3171.636 | 1.426545 | 1.018359 | not applicable |
| rho_incumbent_16 | 2540393 | 3171.664 | 1.426558 | 1.018368 | not applicable |
| rho_incumbent_32 | 2540341 | 3171.599 | 1.426528 | 1.018347 | not applicable |
| rho_scaled_1 | 2967677 | 3705.126 | 1.666499 | 1.189653 | not applicable |
| rho_scaled_2 | 2728437 | 3406.436 | 1.532154 | 1.093749 | not applicable |
| rho_scaled_4 | 2494754 | 3114.685 | 1.400929 | 1.000073 | not applicable |
| rho_scaled_8 | 2540487 | 3171.781 | 1.42661 | 1.018405 | not applicable |
| rho_scaled_16 | 2540537 | 3171.845 | 1.426639 | 1.018426 | not applicable |
| rho_scaled_32 | 2540539 | 3171.846 | 1.42664 | 1.018426 | not applicable |
| rho_pairinv_1 | 2967663 | 3705.108 | 1.666491 | 1.189648 | not applicable |
| rho_pairinv_2 | 2728457 | 3406.461 | 1.532165 | 1.093757 | not applicable |
| rho_pairinv_4 | 2494724 | 3114.647 | 1.400912 | 1.000061 | not applicable |
| rho_pairinv_8 | 2540535 | 3171.842 | 1.426638 | 1.018425 | not applicable |
| rho_pairinv_16 | 2540583 | 3171.901 | 1.426664 | 1.018444 | not applicable |
| rho_pairinv_32 | 2540537 | 3171.844 | 1.426638 | 1.018425 | not applicable |

## Durable artifact and reproduction

The archive is committed beside the existing evidence and registered in
[the manifest](../../evidence/manifest.json):

- File: `ic-reference-qualification-20260925.tar.zst`.
- SHA-256: `a61e6b7c8154dfc594365a044dc74346880018a9c57e14d699ab4c55444950ac`.
- Compressed bytes: 40,999,253; extracted file bytes: 405,702,110; files: 122,855.
- Original gzip profile streams are preserved byte-for-byte. Build caches,
  Python caches and transient operation locks are omitted.
- `initial/`: run 36126765057, implementation head d8daee91959058a10f3c2cfde70e4496af62b05d.
- `final/`: run 36127830907, implementation head 3b1112b65f11c6d3678e8bea26c0a22298e2efa1.
- `full/`: run 36127866931, same final head, including the qualification panel.

Each of these three runs retains 39 IC and 13 rho producer pairs, 17 successful
native driver jobs plus one intentional failure, 15 driver tournament pairs,
and the separate 66-pair qualification wiring control. These are correctness/
integration history, not three independent qualification panels.

From the repository root, on Ubuntu 24.04 with Python 3.12 and zstd:

```sh
python3.12 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-reference-qualification-20260925 --out /tmp/ic-reference-evidence
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/reference-qualification/audit_archive.py \
  --bundle /tmp/ic-reference-evidence/ic-reference-qualification \
  --out /tmp/ic-reference-audit.json
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/reference-qualification/export_report.py \
  --round /tmp/ic-reference-evidence/ic-reference-qualification/full/ic-reference-qualification-36127866931-1/tournament \
  --out research/ic_candidate_tournament_20260915/goal_20260924/reference-qualification/RESULTS.json --verify
```

These commands check stored measurements; they do not rerun measured workers.
Next: bind the qualified sources/settings into the improvement protocol, seal
cross-campaign point exclusions and the familywise confirmation rule, then
execute the bounded rounds with diverse complete pipelines. The active goal
remains unfinished until those rounds produce a qualifying winner or an audited
no-winner decision within the declared budget.

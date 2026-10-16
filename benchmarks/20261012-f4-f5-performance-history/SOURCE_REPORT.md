# Separate F4/F5 plot groups and PDFs

Four groups now retain the algorithms and phases separately: shared Boolean matrix steps, shared Boolean solver calls, complete-call allocation traffic, and F5 criterion construction. The first two each have instruction and diagnostic-runtime line plots. The allocation groups each have cumulative requested-byte and allocation-request plots.

The [combined PDF](reports/run-20261012T0009Z/report.pdf), [graph report](reports/run-20261012T0009Z/REPORT.md), [editable source](reports/run-20261012T0009Z/report.typ) and [native plot data](reports/run-20261012T0009Z/plot-data.json) are retained together. Group PDFs:

- [Boolean matrix steps](reports/run-20261012T0009Z/boolean-matrix.pdf)
- [Boolean solver calls](reports/run-20261012T0009Z/boolean-solver.pdf)
- [Complete-call allocations](reports/run-20261012T0009Z/complete-allocations.pdf)
- [F5 criterion allocations](reports/run-20261012T0009Z/criterion-allocations.pdf)

## Recorded values and scope

Seven existing main-revision snapshots were imported, spanning UTC commit dates 2026-10-08 through 2026-10-11, with revisions `93ff0c72`, `6f81260f`, `8a41eca7`, `0cfba411`, `f88e5909`, `fb926c02` and `e478b952`. Each has one recorded median from five native samples per kernel, plus historical profiler Ir. The snapshots were Python-orchestrated; the exact files, original Git blob IDs and byte sizes are preserved in [the import manifest](records/import-manifest.json), and the graph data retains native SHA-256 hashes, hosts, compilers, harnesses and sample counts.

There is no stored timing spread or exclusive-host eligibility. Historical CI runtime is diagnostic. Historical profiler scope misses some worker-thread work and is not a complete solver operation count. Lines stop at a changed fingerprint, host, compiler, harness or profiler; absent values remain absent. Dates label equally spaced measured revisions, not elapsed time. No new performance measurement, fitted trend or runtime-speedup claim is made.

Allocation history uses the existing [accepted native comparison](../f4_f5_row_build_20261011/records/comparison.json), unoptimized production revision `e957e48c` and optimized revision `2e3dffd`. It has two measured production revisions, eight frozen systems and 24 calls per arm and revision. A/A and the five replication pairs remain controls rather than extra timeline points. Cumulative requested bytes are not peak memory. These workloads are separate from the historical n=17/n=23 kernels; branch allocation results are not inserted as measurements of main.

## Automatic regeneration

The [group registry](../../docs/perf/groebner-groups.json) selects the shared matrix-step and solver cases. The [allocation index](../../docs/perf/f4-f5-allocation-history.json) retains accepted studies and can be appended for later matching native audits. The [Rust tool](../../tools/groebner-history/README.md) prepares plot data without running a benchmark. It revalidates raw allocation records, output/accounting equality, shared inputs, per-call counts, totals and the workload/repetition contract before generating a point. It also calculates provenance hashes and preserves history discontinuities.

The [render-only script](../../scripts/perf/render_groebner_history.py) lays out these prepared values, writes Markdown and editable Typst, and compiles the PDFs. Rendering dependencies are pinned. The [performance-history workflow](../../.github/workflows/perf-history.yml) regenerates every form together, appends links beside the aggregate history, commits the outputs to `perf-history` and uploads them as one artifact. Its PR path validates and renders the new reports. The live workflow runs after publication; local validation is recorded separately and is not a claim that GitHub Actions has passed.

## Validation and delivery

The native package passed **19 tests**, its release build and Clippy with warnings denied. All seven imports matched their original Git blob IDs and byte sizes. Native preparation revalidated the complete allocation records and totals. Workflow YAML parses and both report-generation paths are present. All eight SVGs, five PDFs, sources, Markdown and JSON reproduced byte for byte on a separate render.

Every PDF page was opened and inspected: four single-page group PDFs and all four pages of the combined PDF. Axes, captions, revision labels, condition breaks, tables and source pointers are legible, with no clipping or overflow. PDF inspection and reproduction evidence are retained in `records/`.

The repository library check ran with a 20-minute budget and exited **124** on timeout. Across it and focused reruns, **3,904 distinct library tests passed and 110 were ignored**; one unchanged isogeny test remains unfinished: `cryptanalysis::jv_isogeny_walk::tests::end_to_end_from_a_non_weak_curve_recovers_and_verifies`. The initial parallel run recorded nine socket-fixture failures and two cover-test failures. All socket fixtures passed in a network-enabled rerun (37 passed); both cover failures passed on a serial rerun (14 passed, one ignored). No library source was changed. The raw initial failures, reruns and exact named coverage are preserved. This is **not a full-suite success**. The earlier allocation study retains its original verification record.

Initial tool attempts are retained too: the shared hash tests needed their `hex` development dependency; allocation validation was corrected to use each arm's actual output/accounting schema; and the renderer's Typst table syntax was corrected before successful compilation. The missing Clippy component was installed, and its first formatting warning was fixed. No failed attempt generated an accepted result or changed a frozen measurement.

| Requirement | Status |
| --- | --- |
| Separate F4/F5 matrix and solver groups | Verified complete: separate instruction and diagnostic-runtime series on the shared workloads |
| Complete-call and criterion allocation groups | Verified complete: separate groups, raw outputs and all totals revalidated natively |
| Group and combined PDFs | Verified complete: five PDFs, every page inspected, eight SVGs and editable sources retained |
| Continued regeneration | Implemented and locally validated; GitHub execution occurs after publication and is not claimed as passed |
| Preserve historical measurements and provenance | Verified complete: seven original raw snapshots, Git blob identities, byte sizes and native SHA-256 hashes |
| Repository library verification | Partial: 3,904 passed, 110 ignored, one timed-out unchanged test; initial failures and successful reruns retained |
| PR delivery | Committed update for #1677; publication is the final authoring step |

Canonical graphs checked: the new group outputs and the performance-index documentation change. `docs/index-calculus-scoreboard.html`, `docs/ic/progress-timeline.json`, curve diagrams and the aggregate index remain unchanged: this rendering adds no IC/rho ratio, curve fact, exponent, newly priced phase or aggregate-index value.

Conductor is unavailable; [initial scope attempt](records/conductor.log) and the rendering dependency scope expansion are retained. The requested requirements are frozen in [REQUIREMENTS.md](REQUIREMENTS.md).

## Planned follow-up

The [next-steps roadmap](NEXT_STEPS.md), [editable workflow diagram](next-steps.svg) and [one-page PDF](next-steps.pdf) propose controlled native runtime and peak-memory measurements, worker-inclusive profiling and held-out scaling before choosing another optimization. The unfinished isogeny test remains an independent validation task. These are planned follow-ons; this documentation update launches no new experiment and changes no quantitative graph values. [Planning validation](records/next-steps-validation.md) retains PDF inspection and byte-reproduction checks; the existing library-validation gap remains unchanged.

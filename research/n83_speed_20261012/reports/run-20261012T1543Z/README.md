# Exact scalar and expanded rho follow-up report

**Completed native campaign; descriptive source-linked report.** The actual
first timed child started in UTC minute `20261012T1543Z`. All 139 rho processes
complete and independently verify; no process timeouts, total-step failures,
invalid answers or contended stage
reservations occurred. The final PDF and page-level QA are linked below.

This new report answers the request to make the retained degree-83 rho solver
faster, while keeping the previously tuned Fold/reference-scalar 64-lane DP6
control. It adds a new source-linked report to the comprehensive edition; the
prior 60-page PDF stays read-only and its text and rasterized pages must match
after merging. The old IC/rho results keep their original matched observations.

The only executable Python in this directory lays out native recorded values,
renders documents, or performs document-layout QA. It does not implement a
solver, generate workloads, select a policy, compute a statistic, or replay an
answer. Required native inputs are listed in `native_report_data_contract.json`.

The report has a cover, contents, and eight evidence sections: outcome;
curve and sampling; exact arithmetic; A/A plus arithmetic selection; complete
policy grid; all held-out cases and native intervals; accounting and host;
validation, requirement status and linked reproduction evidence. A failure or
inconclusive comparison remains in the report.

The selected configuration retains reference arithmetic and uses 256 lockstep
lanes, DP6, Fold field, escape off. Native held-out current/selected process
geometric ratio is **1.069651413981749**, target-block 95% interval
**[1.0431783781923283, 1.0892379529034781]**. Mean process costs are
**2.987007664533333 s** versus **2.7805452334000003 s**. All five target means
are lower. The lower interval edge overlaps the descriptive A/A range; this
remains a conditional one-host, five-target result.

Exact Barrett is tested and retained as an unselected candidate. Its nine-case
training mean is 2.923128010666667 s versus reference 2.9030551186666664 s.
All nine exported-semantic pairs match. Internal capped-walk restarts are
retained and their costs remain inside the measured solver intervals. The
report does not claim a scalar
kernel gain or a complete walk trace. Four Lean mathematical declarations
pass their checks; the Rust machine-word bridges remain written/native.

The authoritative native inputs are copied byte-for-byte into `data/`; their
sources and hashes are in `native_report_data_contract.json` and
`report_artifacts.sha256`. No statistic is estimated by the renderer. Rebuild
after timing has ended:

```sh
python3 figures.py
python3 build_report.py
```

The deliverables are [the new addendum](n83_faster_rho_evaluation.pdf),
[editable report source](report_source.md), [six SVG/PNG visuals](figures/),
[layout receipt](render_receipt.json), and [visual review](visual_review.json).
The user-facing combined PDF is
`/workspace/report_output/unified_curve_basis_report_faster_rho.pdf`.

`build_report.py` checks
page bounds and fonts, writes page previews and a contact sheet, and checks
the old 60 pages at text and 72 dpi pixel levels. It does not replace visual
inspection. Commit the source, new PDFs, native input copies, SVGs, 150 dpi
PNGs, layout receipts and visual-review receipt. Identical browsing copies
live in `benchmarks/20261012-n83-exact-scalar-rho/`; the study remains the
source of truth. The first layout attempt's invalid diagram arrow and the
subsequent orphan-source-note pages are retained in rendering attempt logs;
only document layout changed.

The companion [independent scalar audit PDF](../../audit/independent_scalar_audit.pdf)
contains the full written proof and formalization map. The
[fresh-generation report](../../verification/GENERATION_REPORT.md) and
[generation SVG](../../verification/generation.svg) retain the native
exclusion proof. Their source diagrams and the audit PDF are also copied to
the benchmark browsing entry.

Conductor was attempted before edits and is unavailable (`command not
found`). The research-visuals and report-evidence skills govern the evidence
and document workflow. Canonical operation-boundary graphs remain unchanged
unless their admission gates are met; a host runtime improvement does not
change the generic rho boundary or a calibrated operation-equivalent S row.
The checked canonical graphs are `docs/index-calculus-scoreboard.html`,
`docs/ic/progress-timeline.json`, `docs/ic/BOUNDARY_TARGETS.md` and the original
study-local charts. They remain unchanged because this is a rho engineering
comparison, with no new matched IC/rho ratio, calibrated S or admission gates.

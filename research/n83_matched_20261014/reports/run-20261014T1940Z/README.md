# Matched campaign report — 14–15 October 2026

The campaign started at **2026-10-14 19:40:24 UTC**. All 32 processes finished
and independent native replay passes for 218 scalars, 104 IC four-point
witnesses and five same-work A/A pairs. This directory presents native results;
Python performs document layout, plotting supplied values and PDF inspection.
It does not recompute statistics, phase accounting, counters or correctness.

The [17-page addendum](n83_matched_ic_rho_evaluation.pdf) precedes the frozen
earlier 70 pages in the [87-page comprehensive report](unified_curve_basis_report_matched_ic_rho.pdf).
Both contain source links pinned to native-evidence commit
`fb26ccaef5c8f5e97aa81b3503d42e5501b030e4`. The measured source freeze is
`c6e8daafb5a59d27e1be88a2f61f0bd6ee83ef03`; prospective executable/matrix
metadata was published at `ba94c10811ab48da60c9af1824c1eb767a9e17ec` before
the first timed solver child.

All ratios are **IC/rho; below one favors IC**. Actual fresh singleton
whole-process means are IC 1159.349882 s and rho 3.071649 s; the native
eight-pair ratio is 377.435691, with paired 95% interval
[311.812847, 499.010210]. Actual genuine 32-target process means are
1232.975217 s and 78.957400 s; the three-cohort ratio is 15.615702, interval
[14.743543, 17.614848]. These are retained-base, one-guest measurements.

Secondary query means give IC/rho 0.716385, target-block interval
[0.460389, 1.013127]; crossed target/cohort sensitivity is
[0.452168, 1.046078]. Their IC/rho validation/input/output boundaries differ.
The native additive model gives N=1,685, outside executed N=32, with 814 of
20,000 draws nonfinite (fraction 0.0407); this is no observed crossover.
Projected costs may exceed the current 2100-second process cap; larger target
counts require a new prospective equal-budget validation protocol.
Original K600 selection/construction cost remains unavailable. Calibrated S,
ecbench L2, independent-host reproduction and full IC rank/log transcript
validation remain unset. All processes completed under 12 GiB RLIMIT_AS;
the null `memory_limit_hit` indicator does not establish zero memory-cap hits.

All raw attempts, the failed frozen V2 audit and failed original collect
operation remain preserved. The versioned V3 numeric-copy correction and
separate token-preserving collect-only extractor change no frozen measured
bytes and launch no new solver. Native supplement values retain 4,488 charged
internal capped-walk restarts in the 104 main rho answers and 558 in the 10
A/A answers; they are distinct from whole-process failures.

`data/` holds byte-identical copies of authoritative native presentation
inputs, bound by [input_copy_receipt.json](input_copy_receipt.json) and
[native_report_data_contract.json](native_report_data_contract.json).
[report_content.json](report_content.json), [report_source.md](report_source.md)
and [figure_data.json](figure_data.json) are editable presentation source.
The seven figures are SVG plus 150 dpi PNG under [figures/](figures/).

The [render receipt](render_receipt.json) reports zero bounds/glyph issues,
71 addendum links and exact prior-70 page text/pixel preservation at 72 dpi.
Independent final PDF inspection [passes](independent_qa/REVIEW.md) for all
17 new pages. The [QA receipt](independent_qa/qa.json) confirms prior-70 text
and 144-dpi pixel preservation, plus new-17 text/pixel equivalence to the
comprehensive prefix at both 72 and 144 dpi. Earlier layout attempts remain
as layout-only logs/receipts; native results were never recomputed.

To reproduce presentation with the installed ReportLab/Fitz/matplotlib/fonts,
from the repository root:

```sh
python3 research/n83_matched_20261014/reports/run-20261014T1940Z/bind_sources.py
python3 research/n83_matched_20261014/reports/run-20261014T1940Z/bind_native_display.py --statistics research/n83_matched_20261014/reports/run-20261014T1940Z/data/statistics.json --analysis-input research/n83_matched_20261014/reports/run-20261014T1940Z/data/statistics_input.json --replay research/n83_matched_20261014/reports/run-20261014T1940Z/data/independent_replay_v3_receipt.json --supplement research/n83_matched_20261014/reports/run-20261014T1940Z/data/presentation_supplement.json --presentation-metadata research/n83_matched_20261014/reports/run-20261014T1940Z/presentation_metadata.json
python3 research/n83_matched_20261014/reports/run-20261014T1940Z/plot_figures.py
python3 research/n83_matched_20261014/reports/run-20261014T1940Z/build_report.py --prior-report research/n83_speed_20261012/reports/run-20261012T1543Z/unified_curve_basis_report_faster_rho.pdf
```

The default renderer also writes the requested user-facing PDFs under
`/workspace/report_output`; the comprehensive output is copied here for
publication. Final benchmark copies and the complete study archive are owned
by the root workflow, with no new statistical computation in packaging.

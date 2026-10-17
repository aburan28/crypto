# Completed optimized-rho evaluation report

This report renders the completed native campaign: 18 kernel-training solves, 108 policy-training solves, 30 held-out rho solves and five prepared IC target solves. Independent native replay verifies all 161 answers and five IC witnesses, with zero invalid rows. Native selection freezes `fold-L64-DP6-E0` before held-out evaluation; the search is bounded and does not claim a global optimum.

The held-out default/selected fresh-process geometric ratio is **1.756671**, with 95% target-block interval **[1.540157, 2.033353]**. Prepared IC/selected-rho online gives **0.431011 [0.137988, 1.776388]**. All five IC cases are shown, including the 13,034.487993 ms tail. Its arithmetic mean-time ratio **1.542600** points in the opposite direction. There is no definitive overall online winner.

The IC cohort shares one measured preparation across five targets, starting from the retained K600 base; it is not five cold single-target replicas. Original base construction and OS-cache flushing are not measured. The historical five-pair cold ratio is unchanged, and new rho times are never paired with old IC times. Five target blocks follow the declared orbit-exclusion sampling law; three seeds per rho arm do not create fifteen independent targets.

## Artifacts and evidence

- [Addendum PDF](n83_optimized_rho_evaluation.pdf), [editable Markdown](report_source.md), [report content](report_content.json) and [template](report_template.json).
- [Rendering receipt](render_receipt.json), [page contact sheet](preview/contact_sheet.jpg), and [layout QA](preview/layout_qa.json). Combined-document QA is in `combined_preview/`.
- Native recorded inputs are copied without analysis into `data/`; authoritative records remain [in the study](../../records/panel/). The statistics SHA-256 is `6250c42f5a7f6644db5c7c39bd5c712155b2667d081309b5e99202f64f8918f7`.
- `figures/` contains six editable SVGs and rendered PNGs. Quantitative plots retain recorded native values and intervals; no Python estimator, selector, bootstrap or cryptographic calculation is used.
- [Native data contract](native_report_data_contract.json) documents a proposed rendering interface. It is a contract, not the authoritative native schema or a result record.

The PDF distinguishes measured source commits `4c9da3bd` / `4ae7caa2` and immutable executable hashes from the later merged delivery tree. Captured pre-measurement checks remain historical; separate delivery validation belongs in the study's check logs.

## Rebuild

After native timing has ended:

```bash
python figures.py
python build_report.py
```

ReportLab formats source-linked pages; Inkscape rasterizes editable SVG diagrams; PyMuPDF checks page bounds, glyphs and link counts and writes page previews. Every addendum page and all six figures are visually inspected at delivery. The merge reads the previous 48-page report without changing it and checks each historical page's text is preserved.

Default user artifacts are `/workspace/report_output/n83_optimized_rho_evaluation.pdf` and `/workspace/report_output/unified_curve_basis_report_optimized_rho.pdf`. The combined PDF leads with this addendum, followed by the historically frozen prior report. Override the read-only prior source with `--prior-report PATH`; `--addendum-only` skips merging.

Separate post-merge delivery checks pass: release library 3992 passed / 0 failed / 112 ignored in 402.32 s, and the configured example 4 passed / 0 failed in 12.19 s. These checks do not replace the pre-measurement results or measured-source identities.

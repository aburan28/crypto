# Optimized-rho addendum — pending source scaffold

No optimization measurements are supplied or claimed by this scaffold. Do not render or run heavy plotting while a native timed campaign is active.

The intended report is a cover, contents and six sections, approximately 5–8 pages where content permits:

1. New tuning question and the historically frozen default-reference result.
2. Native training grid → frozen selection → held-out evaluation.
3. Exact-field backend and walk-policy fairness, with charged conversion/setup.
4. Native held-out values/intervals and all unsuccessful/capped rows.
5. Unequal stop envelopes and prepared-online versus fresh-process costs.
6. Native replay/validation, evidence register and updated interpretation.

The proposed [native data contract](native_report_data_contract.json) is a rendering contract, not a result or existing native schema. It records the fields needed for the final report. Fill [report_content.json](report_content.json) from [the template](report_template.json) only when native outputs and complete dispositions arrive. `PENDING:` fields block rendering. Keep native estimators, units, confidence scope and intervals verbatim; no Python result analysis is allowed.

Known new stop rules: rho arms share a 60-second process cap; IC has 1800 seconds for preparation. IC's third CLI argument 29 is a rank seed, not a probe budget. Its extraction traverses the full index.states once, with no 20-million-probe cap. Keep all censored/exhausted rows; an eligible online comparison uses uncensored verified complete intervals.

Editable diagrams are SVG source only until timing ends. Quantitative plots are added from recorded native values after the actual result schema is available; no training-best result is substituted for holdout.

After the root agent confirms timing is complete and all pending data is resolved:

```bash
python build_report.py
```

The local ReportLab engine writes the source-linked addendum PDF/Markdown and page/figure QA. The merge opens the prior 48-page PDF read-only and checks that its page text is preserved. It prepends the new addendum and keeps historical bookmarks under a labeled section. It never writes either prior 14-page or 48-page PDF.

Default user artifacts:

- `/workspace/report_output/n83_optimized_rho_evaluation.pdf`
- `/workspace/report_output/unified_curve_basis_report_optimized_rho.pdf`

The prior combined input defaults to `/workspace/report_output/unified_curve_basis_report_updated_n83.pdf`; override with `--prior-report PATH`. `--addendum-only` omits merging. The root agent will assign the final report/run identity once the native campaign's UTC start and durable records are available.

Final closeout requires visual inspection of every addendum page and graph, combined-document page integrity, clean layout/glyph QA and source-linked native evidence. This scaffold has not been rendered or tested against fabricated results.

# Planning artifact validation

Scope: `NEXT_STEPS.md`, its editable SVG, Typst/PDF, the study report link and planning records. No library, benchmark, renderer or workflow implementation changed. No new performance run was launched.

- Typst 0.15.0 compiled `next-steps.typ` successfully. `pdfinfo` reports one A4 page.
- The sole PDF page was rasterized with `pdftoppm` and opened with `view_image`. Initial diagram text overflow was corrected by shortening labels and reducing body font size; the final page was opened again and checked. Text, arrows, source links and the planned-status labels are legible without clipping.
- An independent compilation to a scratch PDF matched the committed PDF byte for byte (`cmp`, exit 0).
- Local Markdown destinations in `NEXT_STEPS.md` and `REPORT.md` exist. The noise-gate link uses the existing `The formula` heading; the baseline/isolation link uses the existing section-10 heading.
- `git diff --check` passed. Documentation and visual changes were reviewed together. Existing quantitative records and graphs were not changed.
- Library tests were not rerun for this documentation update. The prior source validation and its unfinished isogeny-test gap remain in `validation.json`; this is not a new full-suite pass.

PDF SHA-256: `c09a1753fdf9d9af9caebcef48c87e3ba15434007a7f3eb94966b4fdf1a2c3e1`.

Rebuild command, from the repository root:

```sh
/tmp/f4-f5-report-venv/bin/python -c 'import typst; typst.compile("research/f4_f5_history_20261012/next-steps.typ", output="research/f4_f5_history_20261012/next-steps.pdf")'
```

The virtualenv is environment-local; the rendering dependency is already pinned in `scripts/perf/groebner-render-requirements.txt`.

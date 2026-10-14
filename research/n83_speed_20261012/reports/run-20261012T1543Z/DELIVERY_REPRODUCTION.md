# Durable PDF inputs and output

The inspected comprehensive70-page edition is committed here as
unified_curve_basis_report_faster_rho.pdf. The immutable earlier60-page input
is copied byte-for-byte to input_pdf/prior_60_pages.pdf. Their checksums are in
additional_delivery_files.sha256; these additions do not change the original
46-entry report artifact manifest or any PDF bytes.

To reproduce in a separate report-directory copy, use the pinned historical
input explicitly and direct outputs to new paths:

```sh
python3 figures.py
python3 build_report.py --prior-report input_pdf/prior_60_pages.pdf --output /tmp/reproduced_addendum.pdf --user-addendum /tmp/reproduced_user_addendum.pdf --combined-output /tmp/reproduced_comprehensive.pdf
```

Rendering recomputes document layout only. It never recomputes native research
results. A rebuilt PDF can differ in timestamp metadata; original hashes and
all-page visual/pixel validation apply to the delivered PDFs. Preserve earlier
receipts and run reproduction in a separate directory rather than overwrite
the frozen report's metadata or previews.

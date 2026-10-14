# Independent report inspection

This check only inspects document layout, reported native values, citations,
and preservation of the historical PDF. It does not calculate research
statistics. The Python inspection source is limited to PDF extraction,
rendering, hashing and comparison.

The first inspected version has ten new pages and a seventy-page combined
edition. Every new page was read visually at 144 dpi. All sixty appended
historical pages match their original page text, dimensions and complete
144 dpi raster pixels. The automated inspection reports no page-bounds or
replacement/null-glyph problems, and all expected rounded headline values
are present. Native source records were read independently for the reported
ratio, interval, process means, A/A range and all-answer count.

The review requested these changes before final acceptance:

- After multiplying the reciprocal bound by nonnegative `x`, the lower
  inequality must be non-strict because the canonical domain includes `x=0`.
- Current and selected holdout policies use the same reference-scalar
  monomorphization; the two distinct scalar candidates are reference and
  Barrett. The text and figure caption must make that distinction precise.
- The holdout graph legend should move outside the axes because it overlays
  one of Q2's seed observations.
- The campaign-flow diagram's middle title should be shortened to stay
  comfortably inside its box.

The final ten-page addendum was inspected again after all requested corrections
and the native A/B median/minimum supplement. Every final new page was read
visually. All previously noted issues are resolved. The holdout observations
are fully visible, the zero-product proof case is included, the scalar-code
identity is precise, and internal capped-walk restarts are explicitly retained
and charged. The pooled median reversal is displayed and described, alongside
the primary target-block estimate.

Final inspection passes: ten addendum pages, seventy combined pages, no
page-bounds or replacement/null-glyph issues, and all sixty historical pages
again text and pixel identical at 144 dpi using the actual ten-page offset.
The addendum has 56 clickable links, including the native descriptive supplement
and machine-proof status. The headline values, process means, A/A range,
descriptive medians/minima, test counts and all-answer replay counts agree with
their native records. No research values were recomputed by this review.

Accepted addendum SHA-256:
`ceed3b30b7ee40c68d3854f865f72fa9dabb28b73f73437c421673cbb93796b4`.

Accepted combined PDF SHA-256:
`856edd55dd4172374f0f60b23e8605230ce7e63ee41d025c2bc7e9e6d0a8fb0c`.

The pre-polish file hashes, complete link records and pixel hashes are retained
in `qa_before_final_polish.json`. `qa.json` and `inspection.log` bind final
acceptance to the hashes above.

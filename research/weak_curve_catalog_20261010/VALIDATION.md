# Catalogue validation

The catalogue contains 21 source-linked family/condition records, 814 exact
model records and a separate complete 294-class ordinary p7 oracle. Native
generation and verification pass; the local browser and nine-page report have
been reviewed. The family list is extensible and its global completeness status
is `open`.

| Deliverable | Status and evidence |
| --- | --- |
| Family taxonomy | Verified schema, unique IDs, typed relation hypotheses, nonempty properties and local source links in `catalog.rs`; 13 source records and 21 entries. |
| Model index | Verified five immutable SHA-256 input bindings, exact ICV1/model identity, coefficient/modulus shape, field cardinality, trace/order equation, Hasse bound and distinct identities; 814 models. |
| Complete class table | Exact join to the frozen ordinary p7 oracle: 294 traces, 126 cubic-support classes, 210 quadratic-support classes, 88 overlaps, 248 in the union and 46 union zeros. |
| Direct model labels | Complete ordinary p7 invariant images recognize 14 cubic and 10 quadratic records among 277 route/conversion/control models. The sets are disjoint at the model level; class support can overlap. |
| Larger-degree controls | 25 distinct imported models retain their geometric class-witness provenance. Direct membership of these imported models is unclassified in this index. |
| Independent large sources | All 512 original 192–252-bit sources retained, with bounded search receipts and unresolved class labels. The cubic necessary-condition exclusions remain cubic only. |
| Native tests | 37 focused tests pass: 21 library tests and 16 bin tests, including the three catalogue mutation controls. Raw outputs: `native_tests.stdout` and `native_tests.stderr`. |
| Deterministic data and browser | `weak_curve_catalog check .` reconstructs the index, checks the summary, and requires HTML to match the canonical data and template exactly. `check.stdout` records success. |
| Browser | Search, covering-family filter, 512-source filter, model detail, pagination, all 294 classes, 46-zero and 88-overlap filters pass. Screenshots cover desktop, 390/320px mobile and dark theme. No page overflow, script exceptions or external page requests. Source links remain available with JavaScript disabled. Receipt: `browser_review/browser-check.json`. |
| Diagram | Full native WebKit viewport rendered at 1400 by 860 logical pixels; text, arrows and box boundaries visually checked. |
| PDF | Installed Tectonic compiles with empty stderr. All final pages reviewed: pages 1–8 have identical raster bytes to the individually inspected prior render, and final page 9 was inspected after compact clickable citations were added. Nine final page PNGs are retained under `pdf_review/`. |
| Evidence seals | Native `iso1_two_branch_seal` seals the canonical catalogue and this study separately; CI checks both manifests. The parent two-branch evidence remains immutable. |

The first two native compile attempts exposed string escaping and borrowed
boolean mistakes; a later presentation edit exposed a literal-brace format
error. Corrected builds and tests pass. Their stderr/stdout attempts are
retained. Earlier PDF attempts exposed URL and paragraph overflow; the final
compilation has no warnings. The first two LaTeX sources and compile logs
retain those layout diagnostics.

Type I/II and hyperelliptic-locus recognizers remain unimplemented locally.
Binary-descent and order/embedding-degree fixtures with specified subgroups
remain next catalogue additions. Membership in other families, unspecified
subgroups, endomorphism-order conductors and comparative costs remain
unclassified, null or unmeasured as appropriate. The 814 records have distinct
inventory laws and are not a family-density estimate.

This round adds classification and validated evidence joins. It preserves the
numerical IC ratio graphs and prior two-branch census figures because it adds
no new paired timing measurement. Compiler and base-revision metadata are in
`RUNTIME.txt`; source bytes and rendered artifacts are bound by the seals.

Publication extends [PR 1588](https://github.com/aburan28/crypto/pull/1588).
The prior head's full-library execution and broader CI failures are retained
in [publication validation](../iso1_weak_classes_20261007/publication_validation_20261010/README.md).
These remain repository merge gates; current-head checks are tracked on the
PR. Conductor conflict checks allowed the scope; reattaching task T-896
returned server error 500, so its task status is not used as validation.

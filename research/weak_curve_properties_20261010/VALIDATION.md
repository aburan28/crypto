# Property catalogue validation

Validated 11 October 2026. Native arithmetic and independent compiled PARI
agree on all 814 frozen exact models. The resulting sidecar certifies 360
ordinary 2-volcano levels and retains bounded local levels for 454 models.

| Check | Evidence | Result |
| --- | --- | --- |
| Frozen input identity | `models.json` SHA-256 `9738f6bd6c016732f88e4e7979ac09e9aa6fd16eba96220ec8781c0017b6055e`; native input guard | 814 exact models; altered input rejected |
| Native tests | `native_tests_final.stdout`, `native_tests_final.stderr` | 44 tests pass, including seven new field, conductor, twist, torsion and identity controls |
| Portable identities | `identity_certificates.json`, `identity_exact.json` | 13 certificates; 64 saved and 128 fresh evaluations each; exact integer-polynomial expansion agrees |
| Independent arithmetic | `audit_run2/receipt.json`, `stdout.txt`, `stderr.txt`, archived `source/` | 814 unique rows; complete marker; empty stderr; all native properties and supplied root witnesses agree |
| Property regeneration | `check_final.stdout`, `check_final.stderr` | Rebuilt property records, aggregate counts, identity certificates and generated HTML agree |
| Earlier catalogue | Native `weak_curve_catalog check .` | 21 families and 814 models agree; all 512 large-source class decisions retain unresolved labels |
| Earlier evidence | Native seal checks | Earlier catalogue: 55 files; two-branch study: 5,509 files; both seals intact |
| Browser behavior | `browser-review/browser-check.json`, `browser_check.stdout`, `browser_check.stderr` | Filters, pagination, detail fields and exact versus bounded levels pass; no script errors or external requests |
| Browser presentation | Desktop, dark and 390/320-pixel screenshots in `browser-review/` | No page overflow; property table scrolls within its container; source links available with JavaScript disabled |
| Native source formatting | `rustfmt --check --edition 2021` on the four new Rust modules | Pass |
| PDF compilation and visual review | `REPORT.tex`, `REPORT.pdf`, `REPORT.log`, standalone plot and diagram sources | Three pages compile with installed Tectonic; all pages reviewed; formulas, tables, captions and figures readable; no overfull boxes |

The final source files are sealed with their generated data and presentation
artifacts. `EVIDENCE_SHA256.tsv` binds the new study; a separate canonical seal
binds the original family/model indices and the new property JSON and HTML.
The CI workflow rechecks both seals, rebuilds the sidecar in Rust, and independently
reruns PARI, requiring empty stderr and byte-for-byte agreement with audit run 2.

The seven new tests cover finite-field arithmetic and Frobenius, a reducible
modulus refusal, local fundamental-discriminant cases, exceptional-j twist
descent including an obstruction, a changed rational-root witness refusal,
exact floor and torsion bounds, and certificate replay with mutation refusal.
Identity evaluations complement the exact coefficient comparisons. Formal
proof pointers remain null; the torsion/order and cyclic twist propositions
have deductive proofs in the report.

Point counts and traces come from the independently verified parent corpus.
This round verifies structural properties against those fixed values. Prime
characteristics and irreducible moduli are checked in PARI; the native field
checker also uses Rabin irreducibility. Local discriminant information is
computed without factoring its odd part, while PARI independently searches
trace/norm integrality of `(pi-a)/2^k`. Full rational 4-torsion is checked using
root differences natively and the degree-four division polynomial in PARI.

The successful audit's 16,919 ms is a correctness diagnostic on a shared host,
including process launch, under a 300-second cap. Its wrapper and executed
sources are archived. Audit run 1 stopped during stack-setting initialization:
zero accepted rows, `PROCESS_ERROR`, with raw stderr and sources retained.
The wrapper now supplies the stack allocation on the GP command line. No
property row is taken from that failed run.

The frozen input laws differ: 277 p7 route/conversion/control records, 25
constructed higher-degree controls and 512 independently sampled large sources.
All large sources have full rational 2-torsion by their original model law.
Counts describe those inventories. Fields have 17 or 192–252 bits, in degrees
6, 10 and 14. The selected subgroup order, generator and selected-point
containment remain unspecified. Ten p7 ambient orders are 13-smooth; all other
residual cofactors are explicitly unclassified.

The new exact conductor values concern the 2-part of the curve's endomorphism
order. A certified 2-surface model can still have a nonmaximal order at odd
primes. Full endomorphism conductors are null, and the remaining local levels
have explicit bounds. Proper j, proper displayed coefficients and a proper
isomorphic model are independently labelled. Necessary genus-3 filters retain
their scope, alongside the unchanged parent family and class decisions.

This round preserves all earlier IC ratio and census graphs: the new graph
counts exact and bounded structural classifications, with no new paired timing.
Publication status and current-head CI results are recorded in the PR rather
than inserted into this sealed evidence after publication.

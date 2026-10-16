# Validation receipt — 12 October 2026 UTC

The fixed-source importer built `weak_curves.sqlite` and a second full
reimport matched canonical database digest
`e1a5b4d45b680d3e30e856ab0893ef1ef8d218c2238077ca8f2ce4abba178a99`.
It consumed all 19 SHA-256-locked inputs, checked the seven complete exact
domains and three inclusive diagnostic domains, linked 814 models to 814
properties, and required one proof obligation for each of 21 families.
The complete raw stdout, stderr and exit files are in [validation](validation/).
`shasum -a 256 -c ARTIFACTS.sha256` passed for all nine sealed outputs.

| Command and scope | Recorded result |
| --- | --- |
| `cargo test --release --locked --offline` in this study | Exit 0, 3 tests passed; missing exact/diagnostic domain and malformed Boolean controls |
| `cargo fmt --check` in this study | Exit 0 |
| `cargo clippy --release --locked --offline -- -D warnings` in this study | Exit 0 |
| `cargo run --release --locked --offline -- check` in this study | Exit 0; full reimport, table digest, TSV and SVG matched |
| `cargo run --release --locked --offline -- build` in this study | Exit 0; generated database digest matches the checked reimport |
| `lake build WeakCurveDatabase` in `lean/` | Exit 0, 740 jobs |
| `lake env lean AxiomAudit.lean` in `lean/` | Exit 0; all six public declarations list only `propext` and `Quot.sound` |
| `lake env leanchecker --fresh --verbose WeakCurveDatabase` in `lean/` | Exit 0; fresh replay |
| `cargo run --release --offline --example weak_curve_cubic_identity -- check research/weak_curve_database_20261012/cubic_identity.certificate.json` in repository root | Exit 0; `IDC1h3561e7fee541a297` checked at 32 sealed and 32 fresh points |
| `cargo run --release --offline --example weak_curve_cubic_identity -- issue | cmp - research/weak_curve_database_20261012/cubic_identity.certificate.json` in repository root | Exit 0; sealed certificate reproduced byte for byte |
| `tectonic --only-cached --keep-logs REPORT.tex` in this study | Exit 0; two-page PDF, no overfull/underfull or LaTeX warnings; both pages rendered and inspected |
| `cargo test --release --lib --offline` in repository root | Exit 1; 3,933 passed, 10 failed, 112 ignored; raw output retained |

Nine root-library failures report `Operation not permitted` while binding
local sockets in `pollard_collab::{cairn,net}` tests under the sandbox.
The remaining N83 relation-query test fails a scalar equality assertion
at `src/cryptanalysis/koblitz_index_calculus.rs:15515`. The branch changes
no root `src/` file. These outcomes remain failures of the wider test run;
they have not been relabeled as passing checks. The original shell receipt
wrapper also attempted to assign zsh's read-only `status` parameter after
the run; `validation/root_lib.exit` records the observed test exit of 1.

The initial importer build rejected an SQL placeholder-count error and a
property-key mismatch; both were corrected before the recorded full
reimport. Initial Lean attempts exposed an import-path mismatch, a proof
normalization issue and a reserved identifier; the retained final module
builds and passes the axiom and fresh-kernel checks. The first local
Tectonic invocation crashed in macOS system configuration inside the
sandbox, and the cached runtime succeeded under the approved host
execution. A Swift SVG raster attempt hit a restricted module-cache path;
its retry exhausted space on the small system volume. The recorded SVG
remains native, while `plot_recorded.py` rendered the PDF graph from the
checked `CLASS_SUMMARY.tsv` without computing research labels.

The database import, graph render and proof checks used local CPU execution
on the development host. No speedup or construction-cost measurement is
reported. The 192–252-bit model rows are identity and structural records;
selected-subgroup advantage remains `unmeasured` for all 537 such rows.

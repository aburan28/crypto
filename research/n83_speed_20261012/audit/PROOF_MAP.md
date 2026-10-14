# Formal proof map

The project uses Lean `leanprover/lean4:v4.34.1`, compiler commit
`5045d0056413266e57c625dcd7c365b10e377c52`, Linux x86_64. It imports
only bundled `Init`; there are no external Lean dependencies. The official
release archive URL and SHA-256 are in `commands.txt` and
`lean_archive.sha256`. No `sorry`, `admit`, custom claim axioms, native
evaluation or kernel bypasses are used.

| Readable result | Lean declaration | Exact hypotheses/scope | Status |
|---|---|---|---|
| Proposition 1, distributive limb expansion | `N83ScalarArithmetic.limb_product_expansion` | All five natural-number variables; polynomial expansion before floor/carry extraction | Checked local build; axiom/fresh-kernel receipts linked below |
| Proposition 1, high128 reconstruction and word overflow | None | Four limbs below `2^64`; exact floor/carry and Rust word operations | Written proof and independent GMP tests; machine bridge remains open |
| Proposition 2, representable reciprocal via `B-1` | None | `B=2^128`, `2<=m<=2^64-1` | Written proof and independent GMP tests; successor-division/word bridge remains open |
| Proposition 3, residual bound | `N83ScalarArithmetic.barrett_residual_bounds` | All natural `B,m,x`; `B>0,m>0,x<B`; actual `mu=B/m`, `q=x*mu/B` | Checked local build; axiom/fresh-kernel receipts linked below |
| Proposition 3, exact one correction | `N83ScalarArithmetic.barrett_one_correction` | Same general hypotheses; conditional subtraction equals actual `x%m` | Checked local build; axiom/fresh-kernel receipts linked below |
| Proposition 3, instantiated Rust arithmetic | None | Canonical residues and exact high-product/reciprocal bridge | Manual instantiation, worker tests and independent native arithmetic tests; compiled Rust equivalence remains open |
| Proposition 4, canonical modular addition | `N83ScalarArithmetic.canonical_add_one_subtraction` | All natural `m,a,b`; `a<m,b<m`; no assumed conclusion | Checked local build; axiom/fresh-kernel receipts linked below |
| Proposition 4, wrapping word implementation | None | `m<2^64`, `a,b<m` | Written two-case word proof and independent native tests; BitVec bridge remains open |

The checked declarations certify the mathematical kernel statements.
They do not certify the whole Rust solver, its compiled executable,
timings, subgroup-order computation, target generator, or benchmark runner.
Those have their own native tests, exact fixture checks, frozen provenance
and independent answer replay.

`lake_build.attempt1.log`/`axioms.attempt1.log` checked an earlier version
importing bundled Std. Its fresh kernel check was intentionally interrupted
after about a minute to remove the unnecessary broad import. The final
version imports only Init; later receipts apply to that version. Earlier
syntax/algebra simplification failures remain as `lean_build.attempt1.log`
and `lean_build.attempt2.log`; the direct third compile passed. Final
receipt names and hashes are in `verification_receipt.json`.

The repository workflow
[`n83-scalar-arithmetic-lean.yml`](../../../.github/workflows/n83-scalar-arithmetic-lean.yml)
repeats the core-only build, fresh kernel replay and four-declaration axiom
audit on relevant pull requests, `main` changes and manual dispatch. It rejects
incomplete/native-evaluated proof sources and unexpected transitive axioms.
The local workflow syntax and shell guards are checked separately; no remote
CI execution or CI outcome is claimed by this audit.

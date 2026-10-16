# Native compact-CM implementation and controls, 2026-10-16 UTC

The requested larger-degree campaign now has an installed `isogeny compact-cm`
entry point for a separate known-endomorphism category on the exact registered
secp256k1 source. This implementation round executed fixed controls only. The
proposed target degree 2^8191-1 remains unexecuted, awaiting a clean executable,
immutable custody and its separately frozen experiment approval. Main expanded-map
challenge scores and its 54 separately ranked source models are unchanged.

| Requirement | State and evidence |
|---|---|
| Dedicated native binary | Implemented in the existing isogeny executable; no Cargo at runtime; CLI integration runs with the temporary directory as working directory |
| Prime degree from unreduced expression | Verified on controls with degrees 7, 31 and 127; exact signed integer coefficients and both exact integer identity checkers |
| Independent replay | Separate repository limb integers/general modular reduction and fixed-limb affine group arithmetic; producer uses num-bigint, Mersenne folding and Jacobian scalar multiplication |
| Machine protection | Public construct/verify use existing native supervisor, 1800/3600-second maxima and 8 GiB sampled RSS; positive RSS acknowledgement precedes mathematical work |
| False-positive controls | LL M11 composite with prime exponent11; composite exponent9; unsupported exponent2; altered certificates, wrong beta/lambda pairing, integer errors invisible modulo the IDC1 field, signed/zero fixtures and degree alias |
| Larger target | Not run; proposed M8191 requires separate frozen scientific protocol |
| New cross-platform release | Not verified; existing ReleaseMe configuration builds/tests/packages the same binary, while PR #1976 carries this implementation |

The full native test suite passed: 3 field-check tests, 81 main-binary tests,
4 registry-model tests, the existing challenge CLI integration and the new
compact CLI integration, 90 total with zero failures. The final run exit code is
0 in native-tests-final-retry1.exit-code.txt. Retained native controls are under
native-controls-final-retry1/. Their original commands, stdout/stderr, resource
receipts, mutated certificates, accepted/rejected records and exact point images
remain intact. tested-source.sha256 binds the tested source files.

Earlier passing controls and native test logs are retained. The failed
native-tests-final.log/exit-code.txt captures a test-declaration syntax error
introduced during editing; the corrected fresh retry passed. No scientific
target was launched by that failed compilation. The required root release-library
gate passed separately at unchanged root source: 4148 passed, zero failed,
117 ignored; see ../root-release-retry2.log. Root src and Cargo files are unchanged.

The degree identity's evaluated supplement is IDC1hf56b3716b0a7eb8b in
../../../docs/identity-certificates/cm-eisenstein-norm.json, with its derivation and
trust basis in the adjacent Markdown file. Its finite-field evaluations are
supplementary; the two native exact integer coefficient checks are separate.
The geometric degree conclusion uses the reviewed dual-norm theorem and the
unreduced expression, never rational-point agreement alone. The degree-alias
control demonstrates that distinction concretely.

The exact source parameters are SEC2 secp256k1 with cofactor1 and source
icv1-fp256-t432420386565659656852420866390673177327-76dadd18. Native validation checks
the tuple against the pinned registry, generator/order and unique Hasse multiple.
Primality of p,n is explicitly trusted from SEC2; there are no self-contained
p,n primality certificates here. Target-degree primality is independently replayed.
Pure Frobenius and 1-Frobenius remain theorem-level controls outside this canonical
Mersenne/Cornacchia executable schema; they are not reported as executed inputs.

Constructor and verifier retain native worker receipts. The input is bounded at
2 MiB and its hash binds the same bytes parsed for verification. CLI evaluate
rechecks the certificate, with separately named geometric and rational-subgroup
semantics; constructor/verifier watchdogs do not extend to that standalone evaluator.
All timings in native receipts are operational observations on a shared Apple
M4 Pro host, with no performance comparison or speedup claim. The diagram replaces
a comparative timing graph because this is correctness validation, not a paired
performance experiment. IC/rho, performance charts and canonical challenge graph
values are unchanged.

Mathematical review and implementation authorization: Crypto Autoresearcher
PR #2221, TASK-20261016-a82f9e and DEC-20261016-3fb077. The independent reviewer’s
served model/effort identifiers were unavailable and recorded as unverified.
The implementation is published through https://github.com/aburan28/crypto/pull/1976.

The one-page report.pdf was rendered with Typst 0.15.1 and visually inspected
in full. The precommit control binary's real dirty-tree provenance and SHA256
are retained in tested-binary-version.json and tested-binary.sha256; clean
scientific execution will use a separately built and pinned committed executable.

#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)
#align(center)[#text(size: 18pt, weight: "bold")[Native compact-CM certificates]]
#align(center)[Implementation and controls · 16 October 2026 UTC]

The installed `isogeny` binary now exposes compact-CM construction, fresh independent verification and named geometric/subgroup evaluation. This round validates controls; the proposed degree 2^8191−1 has not been executed. Compact known-endomorphism self-maps form a separate category from the main expanded-map challenge.

#image("controls.svg", width: 100%)

= What ran

The complete dedicated-tool suite passed 90 tests with zero failures: 3 field checks, 81 main-binary tests, 4 registry-model tests and two CLI integrations. Fixed prime-degree controls are 7, 31 and 127. The installed CLI is invoked from outside the checkout. Its positive construction/replay, malformed certificates, resource receipts and point images are retained under `native-controls-final-retry1/`; the final log and exit-code file bind the validation outcome. Root release-library validation separately passed 4148 tests, with 117 ignored and zero failed, at unchanged root source.

The producer uses num-bigint, Mersenne reduction and Jacobian scalar multiplication. The independent verifier uses repository limb integers with general modular reduction, independent Cornacchia replay, dense exact integer identity expansion, and fixed-limb affine point arithmetic. It checks eight fixed inputs, subgroup action, group compatibility and both dual compositions. The public constructor/verifier use native watchdogs and sampled RSS protection; no performance gain is claimed.

= Proof and rejection boundaries

Evaluated supplement `IDC1hf56b3716b0a7eb8b` is accompanied by exact integer coefficient checks. A deliberately wrong integer coefficient invisible modulo the supplement's field is rejected. Wrong beta/lambda pairing is rejected. Adding the subgroup order to a coefficient preserves rational-point images while changing the geometric norm, so its original degree certificate fails. Prime exponent 11 with composite M11 and composite exponent 9 are distinct negative controls. Earlier compile/test logs are preserved.

The exact SEC2 secp256k1 tuple is checked against the pinned registry; p,n primality is explicitly trusted from SEC2, rather than self-contained here. The new degree's primality has its own independent Lucas–Lehmer replay. The reviewed dual-norm theorem establishes geometric degree from unreduced integers; point agreement cannot do so. Frobenius fixtures outside this canonical schema remain theorem-level controls.

= Provenance and next gate

Source hashes, logs and the editable SVG are beside this report. Mathematical review and implementation handoff are in #link("https://github.com/aburan28/crypto-autoresearcher/pull/2221")[PR2221]; implementation publication is #link("https://github.com/aburan28/crypto/pull/1976")[PR1976]. The exact identity and proof basis are in `docs/identity-certificates/cm-eisenstein-norm.{json,md}`. The next scientific gate requires a clean executable, immutable custody and a frozen target protocol. No worldwide record, new expanded-map value or comparative timing result is established by these controls.

The clean retake uses source `9a2e133099` and executable SHA256 prefix `d7e53928e57b`, with dirty=false. One complete native CLI integration passed (exit 0). Raw receipts are under `clean-build-9a2e133099/`. Remote custody awaits explicit upload approval after automatic review rejected the external transfer; no target launch is authorized yet.

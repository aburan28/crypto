#set page(paper: "us-letter", margin: (x: 0.65in, y: 0.58in))
#set text(font: "DejaVu Sans", size: 9.4pt)
#set par(leading: 0.42em)
#align(right)[#text(size: 8pt, fill: rgb("#52687c"))[TASK-20261016-de577e · 2026-10-16]]
= Fixed-modulus replay: source and engineering controls

The goal is larger independently verified prime-degree isogenies on all 54 standard prime-field boards, ranked separately. This source change reuses one kernel reciprocal across selected torsion and subgroup composition. The candidate stays off by default. It establishes no speed comparison or new admitted degree.

#image("architecture.svg", width: 100%)

#table(columns: (2fr, 1fr, 2.3fr), inset: 5pt, stroke: 0.4pt + rgb("#b5c4d0"),
  [*Gate*], [*Outcome*], [*Frozen evidence*],
  [Native reference A], [128 pass / 0 fail], [implementation-02/native-reference-tests.log],
  [Native candidate B], [128 pass / 0 fail], [implementation-03/native-candidate-tests.log],
  [Root release library], [4148 pass / 0 fail], [implementation-02/root-lib.log; 117 ignored],
  [Affected research harnesses], [12 pass / 0 fail], [implementation-02/research-native-tests.log],
  [Retained reference bodies], [8 byte-identical], [reference-functions.json],
  [Public mathematical JSON], [3 pairs identical], [development-semantic-comparison.json],
  [Adversarial controls], [9 pairs + 2 cache negatives], [all rejection status/reason/score equal],
)

Both native suites have one ignored child helper which a parent test invokes. Tests include public P-192/13 and P-521/7 fallback controls, exact accepted P-224/1471 cache coverage, and recomputed-hash bad-kernel and bad-numerator controls. No mathematical fields were excluded from comparison. All twenty transports remain required by the shared admission path.

Host: Apple M4 Pro, Mac16,7, 14 logical CPUs, 51,539,607,552 bytes RAM, macOS 26.6 build 25G70. This is a shared host. Test and resource times are operational receipts, not isolated performance measurements. Toolchain and host files are retained in implementation-01.

#pagebreak()
= Correctness, provenance and next gates

The context borrows the exact field and trimmed modulus. For length L at least 128 it retains L inverse coefficients, including zeroes; each eligible reduction uses its required prefix. The predicate L at most dividend length less than 2L and the quotient/remainder degree assertion are preserved. Zero or trailing-zero moduli are rejected. Original divrem, rem, mulmod, gcd, invmod, dense/selected division and composition bodies are byte-identical to base1180.

Cached reduction is checked against both old reciprocal reduction and schoolbook division, including quotient identity. Controls cover L=127/128/129/130, operand lengths47/48/49, precision endpoints and doubling boundaries, nonmonic/sparse/repeated-root moduli, interleaved fields/contexts, and a corrupted reciprocal. Recurrence is checked against dense division; composition against old baby/giant composition and Horner evaluation.

*Proof map.* Native certificate IDC1h41e9973f839c85bf is a finite evaluation supplement over p=2^61-1. Lean4.35.0 checks the integer scalar Newton residual identity and the integer-coefficient convolution precision lemma, with only propext and Quot.sound axioms. Full formal-series uniqueness and polynomial quotient-reversal Lean bridges remain open. The written proof and exact scopes are in docs/identity-certificates/fixed-modulus-newton.md. No fresh lean4checker was available.

*Failures retained.* The first native compile lacked aliases needed by included tests; its source snapshots and raw log remain in implementation-01. The first Lean check named an unavailable lemma; its failed source/log are also retained. Both repairs passed fresh checks. A JavaScript report-preparation call failed before execution or writes and was corrected.

*Resource contract.* The new research companion has a distinct fixed-modulus-bench worker shape. Parent acknowledgement requires positive RSS before worker parsing/algebra. Existing 8GiB, 25ms, 100ms-grace, done and deadline behavior remains. Lost/zero RSS, acknowledgement, parse failure, writer failure and timeout controls passed. The future analyzer binds native outputs, resource receipts and OS real/user/system timings; CPU and wall are separately adjudicated.

*Remaining source gate.* The source commit is followed by clean matched A/B builds, installed-byte equality, retained installed controls, immutable custody and independent source/resource review. Their later receipt must identify the precise clean source commit. The binaries used above are development builds with dirty=true and are not scientific campaign executables.

*Execution boundary.* No benchmark campaign, sealed50341 replay, larger construction or compact target ran. Scientific comparison requires a later frozen protocol. Full source plan: PLAN.md. The prior M44497 compact archive and all54 expanded boards are unchanged. No canonical degree or speed graph changes; this diagram records the implementation only.


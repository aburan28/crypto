# Independent scalar and protocol audit

The arithmetic audit supports correctness of the new exact scalar worker
under its canonical-input contract. Correctness alone does not establish a
runtime gain. The campaign's arithmetic selection retained the reference
scalar kernel when the new Barrett kernel did not improve measured cost.
All 139 campaign solves completed and passed independent replay. The complete
policy and holdout records control the measured speed claim.

- [Readable arithmetic proof](ARITHMETIC_PROOF.md) and
  [formalization map](PROOF_MAP.md).
- [Lean source](ScalarArithmetic.lean), [axiom audit](AxiomAudit.lean),
  pinned [toolchain](lean-toolchain) and [core-only build](lakefile.toml).
- [Pre-launch protocol audit](PROTOCOL_AUDIT.md).
- [Five-page audit PDF](independent_scalar_audit.pdf),
  [proof and evidence diagram](arithmetic_proof_map.svg),
  [150 dpi rendering](arithmetic_proof_map.png) and
  [PDF layout review](pdf_layout_qa.json).
- [Independent native scalar receipt](cpp_check.json),
  [limb polynomial certificate](limb_expansion.certificate.json) and
  [certificate replay receipt](certificate_check.json).
- [Final arithmetic/formal verification receipt](verification_receipt.json)
  and [campaign answer replay](../verification/final_receipt.json).

The first missing-Boost-header failure and the Lean simplification failures
are retained. The broader Std fresh-kernel attempt was intentionally
interrupted; the subsequent core-only Init build, axiom audit and fresh
replay passed before any timed campaign runs. The fresh kernel receipt
prints `replaying ScalarArithmetic with --fresh` and returned exit 0.

No audit build, native test, proof check or report rendering runs concurrently
with the timed campaign. Text-only audit and document sources can be prepared
while the benchmark reservation is active.

# First Linux v2 outcome: versioned archive replay

Source: [PR #921](https://github.com/aburan28/crypto/pull/921) exact frozen head `d70396a649cf4a12ece39f01bf60f121c6b2a9ae`, Actions run `36534548089`, artifact `11017624684`. This successor contains the original 3,950,207-byte ZIP without repacking (`SHA-256 b719f83d46057e45f83aadffb83287f9c31168afdf985f7beb3b4c69ed1e1d76`). Its receipt and manifest hashes match `FROZEN.json`; the manifest binds 113 payload files totaling 7,444,986 bytes.

| Gate | Observed decision | Meaning |
|:--|:--|:--|
| Original #921 `ci_replay.py` | `FAIL_MANIFEST_ORDER_ONLY` | It compares producer Path-component order with string order and stops at the inventory check; it did not reach semantic admission. |
| Versioned successor source gate | PASS | The copied verifier differs only in its source-directory binding and `Path.parts` manifest comparator. |
| Archive mutation controls | PASS | Both an omitted path and an invented path are rejected with internally updated manifest/receipt digests. |
| Corrected independent replay | `PASS_INDEPENDENT_REPLAY` | Original static freeze/release checks, toy exhaustive relations, 20 n13 SAT/DRAT queries, and n131 canonical DIMACS clauses/output replayed from the immutable ZIP. |

The corrected n131 CNF SHA-256 is `33efea4b5a6d004f693ccb1e1485881363011b4f640be24da62144ed60b273d2`. The toy replay covered n=2 and n=3, with parent-row digest `a2b760b72e477cb9ab4d8409e71102cd6a2a09db0b8a5e8f5ad170d1a1102da9`. Local Python 3.12 replay passed in 7.2 seconds. Hosted exact-head replay and independent review are separate admission gates tracked on this PR.

This is a verifier correction for the one existing stage artifact. It does not alter #921's frozen runner or perform another measured child. The archived receipt's `PASS` and this successor's semantic `PASS` must be reported separately from the original verifier failure; none establishes a PDP, recovered logarithm, or matched-rho speed result.

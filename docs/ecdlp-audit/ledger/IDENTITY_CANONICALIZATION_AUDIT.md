# WS3 canonical JSON identity-contract audit

**Status:** non-authoritative audit note; proposed decision only

**Observed:** 2026-10-14

**Audit basis:** working tree rooted at
`1e447956d3616092a223264f1c6c31eff8ac9979`; the native planning validator
was an uncommitted file in the same change set when this note was prepared.

This note records a concrete identity hazard found while beginning WS3. It does
not create the authoritative audit ledger, promote the planning mock, or close
`WS3-G1`. The existing
[`authoritative-ledger-identity`](../ECDLP_AUDIT_PLAN.mock.json#L706)
obligation remains open: its exit evidence still requires an authoritative
schema, frozen canonicalization vectors, collision/alias tests, and a reviewed
decision.

## Scope and result

| Requirement | Status | Evidence or remaining gap |
| --- | --- | --- |
| Determine whether repository identity producers agree on canonical JSON | verified for the inspected implementations | They do not: one codec escapes every non-ASCII scalar while the documented EC1 codec and two other Rust codecs emit raw UTF-8. Exact references are below. |
| Pin a proposed WS3 behavior | proposed, not approved | Use sorted, compact, raw-UTF-8 JSON with no floats under a new versioned audit-identity codec. Do not rewrite existing identity namespaces. |
| Define cross-language conformance evidence | specified, not implemented | Required positive and negative vectors are listed below. No frozen WS3 vector artifact or independent replay receipt exists yet. |
| Establish authoritative ledger identity | open | Governance review, a strict ingest parser, a shared implementation, migration receipts, and independent conformance replay are still missing. |

This was a static source audit plus byte-level SHA-256 reproduction. It was not
a ledger generation run, scientific experiment, curve assessment, or
cryptographic performance measurement.

## Observed divergence

The repository contains two incompatible canonical encodings for the same
parsed JSON value.

| Source and identity domain | Observed behavior |
| --- | --- |
| [`docs/curve-identities.md` lines 31–37](../../curve-identities.md#L31-L37) | Defines the cross-repository EC1 convention as sorted-key, compact UTF-8 JSON with `ensure_ascii=False`, no floats, significant list order, and exact preimages. It requires a versioned adapter when an existing producer differs. |
| [`tools/curve_identity.py` lines 11–29](../../../tools/curve_identity.py#L11-L29) and [`research/ic_candidate_tournament_20260915/identity.py` lines 16–34](../../../research/ic_candidate_tournament_20260915/identity.py#L16-L34) | The historical reference implementations emit non-ASCII characters as raw UTF-8. |
| [`src/cryptanalysis/ecbench/claim.rs` lines 45–112](../../../src/cryptanalysis/ecbench/claim.rs#L45-L112) | Implements the raw-UTF-8 EC1/IC1-compatible codec separately and pins `é` as raw UTF-8 in its test at [lines 1940–1947](../../../src/cryptanalysis/ecbench/claim.rs#L1940-L1947). |
| [`src/cryptanalysis/pollard_collab/cairn.rs` lines 79–157](../../../src/cryptanalysis/pollard_collab/cairn.rs#L79-L157) | Uses raw UTF-8, escaping only JSON syntax and C0 controls. Its frozen cases include Latin, CJK, astral characters, control characters, and Unicode keys at [lines 1040–1074](../../../src/cryptanalysis/pollard_collab/cairn.rs#L1040-L1074). Cairn remains its own protocol and namespace; matching string behavior does not make it a WS3 authority. |
| [`src/bin/ecdlp_audit.rs`](../../../src/bin/ecdlp_audit.rs) (`expected_plan_uid`, `canonical_utf8`, and `identity_codec_uses_sorted_raw_utf8`) | The new non-authoritative plan validator deliberately uses raw UTF-8 and binds the complete plan value except the top-level `plan_uid`. Its test covers a non-ASCII key and value. The replaced planning validator at commit `1e447956…`, `docs/ecdlp-audit/verify_plan.py:26-33`, also used `ensure_ascii=False`, but hashed a narrower skeleton. |
| [`src/cryptanalysis/ecbench/canonical.rs`](../../../src/cryptanalysis/ecbench/canonical.rs) | The historical codec escapes every non-ASCII scalar as lowercase `\uXXXX`, using UTF-16 surrogate pairs for astral characters; its test fixes `é` as `\u00e9`. At the pinned base revision its header incorrectly claimed EC1 compatibility. The current worktree corrects that documentation without changing legacy bytes or identities. |

At the pinned base revision, the header of `ecbench/canonical.rs` said that it
implemented the convention in `docs/curve-identities.md`, but those byte rules
disagree for any non-ASCII key or value. The current full-record planning UID
contains the multiplication sign in `inventory_freeze.unit_expression`, so it
now exercises the raw-UTF-8 branch as well as the dedicated Unicode unit test.
That improves this mock's regression coverage but does not make the historical
codecs interchangeable or close the broader vector and ingest obligations.

### Reproducing the mismatch

For the parsed value `{"a":[1,-2],"b":"é"}`:

| Codec | Canonical text | UTF-8 bytes (hex) | SHA-256 |
| --- | --- | --- | --- |
| proposed raw UTF-8 | `{"a":[1,-2],"b":"é"}` | `7b2261223a5b312c2d325d2c2262223a22c3a9227d` | `2e93e78b27ecfd1e166b33fcad5d2770861401c7c10843ce8f941f104d79b69f` |
| current ecbench escaped-Unicode codec | `{"a":[1,-2],"b":"\u00e9"}` | `7b2261223a5b312c2d325d2c2262223a225c7530306539227d` | `42361525cbf563fd2b849c56c136bb870170f87e9721a958724b404c6b680a31` |

Reproduction commands:

```sh
printf '%s' '{"a":[1,-2],"b":"é"}' | sha256sum
printf '%s' '{"a":[1,-2],"b":"\u00e9"}' | sha256sum
```

The hashes differ because SHA-256 receives different bytes. This is an
encoding-domain mismatch, not evidence of a SHA-256 collision.

## Proposed WS3 codec: `ecdlp-audit.cjson-utf8/v1`

The proposed decision is to follow the already documented EC1 behavior and
pin raw UTF-8 for new WS3 identities. It is intentionally a proposal until the
WS3 and WS8 review gates approve it.

1. Parse well-formed UTF-8 JSON before hashing. Hash the canonicalized value,
   never the source spelling or a pretty-printed document.
2. Admit only `null`, booleans, integer tokens, Unicode strings, arrays, and
   objects with string keys. Reject floats, exponent notation, NaN, infinity,
   duplicate object keys, malformed UTF-8, lone UTF-16 surrogates, and negative
   zero at the ingest boundary. Integers must not be rounded by the host
   language; implementations need arbitrary precision or a schema rule that
   represents wide integers as validated decimal strings.
3. Sort object keys lexicographically by Unicode scalar-value sequence. For
   valid Unicode this agrees with lexicographic UTF-8 byte order. Preserve
   array order.
4. Emit no insignificant whitespace, byte-order mark, or trailing newline.
5. Emit non-ASCII Unicode scalars directly as UTF-8. Escape `"`, `\`, and
   U+0000–U+001F; use `\b`, `\t`, `\n`, `\f`, and `\r` for their five short
   forms and lowercase four-digit `\u00xx` for the remaining C0 controls. Do
   not escape `/`, DEL, U+2028, U+2029, BMP non-ASCII, or astral scalars.
6. Perform no Unicode normalization. NFC and NFD strings remain different
   values and therefore normally have different identities.
7. Compute lowercase hexadecimal SHA-256 over the exact canonical UTF-8 bytes.
   The schema must bind the record kind, identity-schema version, and codec
   version in the preimage and in the stored identity envelope. A digest alone
   is not a cross-domain join key.

The authoritative schema still has to decide the exact domain-separated
preimage and URN grammar for unit, assessment, evidence, implementation, and
correction identities. Until then, `ecdlp-audit.cjson-utf8/v1` is a descriptive
proposal, not a registered codec or authority.

### Existing identities are immutable

Adopting the proposed codec must not recode an existing ecbench, Cairn, EC1,
ICV1, IC1, workload, run, or archived evidence identifier. A migration adapter
must retain:

- the original namespace, codec/version, full digest, and exact preimage;
- the separately computed WS3 identity and its codec/version;
- an explicit typed relation explaining whether the records are identical,
  equivalent under a reviewed adapter, or merely associated; and
- a receipt showing both byte strings and digests when the codecs differ.

ASCII-only preimages may produce equal digest bytes under both codecs, but that
does not make the codecs interchangeable. Equal digests in different identity
domains require the domain metadata for interpretation.

## Required cross-language conformance vectors

The reviewed vector set must be a frozen, hash-addressed artifact. Every case
must carry an ID, input UTF-8 bytes in hex, expected accept/reject result,
canonical UTF-8 bytes in hex when accepted, expected lowercase SHA-256, and the
codec/version. Literal JSON alone is insufficient because editors can
normalize or re-escape it.

At minimum, the frozen set must cover:

| Class | Required case and expected result |
| --- | --- |
| Sorted ASCII baseline | `{"b":2,"a":1}` canonicalizes to `{"a":1,"b":2}` with SHA-256 `43258cff783fe7036d8a43033f830adfc60ec037382473548ac742b888292777`. |
| Raw BMP value | `{"a":[1,-2],"b":"é"}` canonicalizes to the raw form and digest `2e93e78b27ecfd1e166b33fcad5d2770861401c7c10843ce8f941f104d79b69f`. It must not produce the escaped-codec digest above. |
| Unicode key ordering | A permutation of `A`, `Z`, `a`, `b`, and `é` must sort by Unicode scalar sequence and emit `é` raw. Include non-ASCII at more than one nesting depth. |
| Astral and CJK values | `{"astral":"😀𝔸","cjk":"日本語"}` must remain raw and hash to `d4ecaf6f3034360d1ec80794c47a6ad0b24397fdf9022f4120003715b5d002e8`; no UTF-16 surrogate escape may appear. |
| Escapes and controls | Cover quote, reverse solidus, solidus, every short control escape, another C0 control, DEL, U+2028, and U+2029. Pin exact lowercase escape bytes. |
| Equivalent source spellings | Raw `é` and source escape `\u00e9`, plus raw astral text and a valid surrogate-pair source escape, must parse to the same value and produce identical raw canonical bytes. |
| No normalization | NFC `{"label":"é"}` hashes to `59a69e5a8f92914b2096dd10b10431fe5c21baa7990508f95718e2f40fc718db`; NFD `{"label":"é"}` hashes to `bf1c0db8e2dcab417156b09650994aeebd8c0275198748bccb42f0b97de669c2`. Both must remain distinct. |
| Integer boundaries | Cover zero, positive/negative values, `i64`/`u64` edges, values above `2^64`, and values above JavaScript's exact-integer range. No implementation may round or switch to exponent notation. |
| Rejected number forms | Reject `-0`, decimal fractions, exponent notation, NaN, and infinities before identity calculation. |
| Rejected structural input | Reject duplicate keys at every depth, malformed UTF-8, lone surrogates, trailing data, and non-JSON input. |
| Stability properties | Reparse-and-recanonicalize is byte-idempotent; object insertion order never changes output; array permutation normally does. |

The authoritative Rust implementation must pass the frozen vectors. At least
one independently written non-Rust implementation must reproduce every
accepted byte string and rejection before WS3 promotion. The repository's
historical Python references may be checked as provenance, but no new Python
validator or execution path is required or proposed. Independent replay must
record implementation revision, toolchain, vector-artifact digest, command,
exit status, and result digest.

## Open obligations

| ID | Status | Exit evidence |
| --- | --- | --- |
| `WS3-ID-01` | open | WS3 and WS8 reviewers approve or replace the codec rules, domain separation, identity envelope, and URN grammar in a versioned decision record. |
| `WS3-ID-02` | open | A shared native canonicalizer and strict ingest parser replace audit-specific duplication; duplicate keys, negative zero, floats, malformed Unicode, and wide integers have fail-closed tests. |
| `WS3-ID-03` | open | The full positive/negative vector artifact is committed, content-addressed, and replayed by the authoritative Rust implementation. |
| `WS3-ID-04` | open | An independent non-Rust implementation reproduces the frozen vectors and emits a signed or hash-addressed replay receipt. |
| `WS3-ID-05` | open | Every repository identity namespace is inventoried with its actual codec, preimage schema, digest format, and Unicode behavior; misleading same-convention claims are corrected. |
| `WS3-ID-06` | open | Migration adapters preserve legacy identifiers and emit dual-byte/digest receipts for divergent codecs; ASCII equality and non-ASCII divergence receive explicit alias/collision tests. |
| `WS3-ID-07` | open | Authoritative unit, assessment, evidence, implementation, and correction schemas bind identity kind, version, full digest, exact preimage, and supersession semantics. |
| `WS3-ID-08` | open | The master-ledger semantic verifier rejects unknown codecs, codec/digest mismatches, digest-only cross-domain joins, duplicate latest revisions, and unreviewed adapters. |

Until all applicable obligations close, the correct status is **identity
contract proposed; authoritative ledger identity not established**.

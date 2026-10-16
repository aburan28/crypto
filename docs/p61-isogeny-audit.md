# External p61 census ledger and uniqueness audit

`p61_isogeny_audit` supplies the global control that the bounded-memory class
producer deliberately omits. It seals an ordered segment ledger from immutable
generation and independent-verification receipts, then extracts every canonical
j-invariant into fixed-width records, performs an external merge sort, and
reports exact duplicates.

The sorted index format is:

```text
offset  size  meaning
0       8     ASCII P61JIDX1
8       8     little-endian record count
16      16*N  repeated (j u64 little-endian, node_index u64 little-endian)
```

Each temporary chunk and the final index is create-only, independently hashed,
and retained. The audit checks the following before a j-invariant enters the
sort:

- stored segment, checkpoint, generation-receipt, and verification-receipt
  digests;
- agreement between generation and independent-verification envelopes;
- source commit and class-number-receipt binding;
- contiguous checkpoint hashes and x ranges;
- exact `(x,y,z) -> node_index` mapping;
- exactly one node for every local segment coordinate;
- the raw JSONL record-chain digest and terminal summary;
- eighteen terminal closures when the ledger reaches the row order.

The class-wide success state is available only for coverage `x=0..19894128`,
exactly 358,094,304 unique j-invariants, all terminal closures, and a ledger
with a validated, dedicated PARI flag-0 certification receipt. Existing
segments remain bound to the immutable provisional preflight receipt; the
certificate separately binds that exact receipt by SHA-256 and upgrades its
class-number status. This dual binding avoids rewriting or relabelling prior
segment checkpoints.
Otherwise a valid prefix is reported as `partial_census_pass`, or complete data
with a conditional class number is reported as
`full_data_pass_certificate_blocked`.

Build a ledger from a release directory:

```sh
target/release/p61_isogeny_audit build-ledger \
  --generation-receipts-dir receipts \
  --verification-receipts-dir verification \
  --class-number-receipt ../p61-class-preflight-v1/receipts/preflight.json \
  --class-certification-state conditional_grh \
  --auditor-source-commit <40-hex-commit> \
  --output ledger.json
```

After `bnfcertify(bnf)` with PARI flag 0 returns `1`, build an unconditional
ledger by supplying both pieces of evidence:

```sh
target/release/p61_isogeny_audit build-ledger \
  --generation-receipts-dir receipts \
  --verification-receipts-dir verification \
  --class-number-receipt ../p61-class-preflight-v1/receipts/preflight.json \
  --class-number-certification-receipt ../p61-class-certification-v1/receipts/certified.json \
  --class-certification-state unconditional \
  --auditor-source-commit <40-hex-commit> \
  --output ledger-certified.json
```

The auditor rejects an `unconditional` declaration without the dedicated
receipt. It also checks the frozen field, trace, discriminant, class number,
cyclic factors, all three class-number computations, `bnfcertify` flag and
result, source identities, and the provisional-receipt digest. A certificate
cannot be substituted for the receipt hash already embedded in the segments.

Audit it with bounded RAM. Five million records require about 80 MB for the
pair buffer, plus sorting overhead:

```sh
target/release/p61_isogeny_audit audit-ledger \
  --ledger ledger.json \
  --work-dir audit/chunks \
  --sorted-index-output audit/p61-j-index.bin \
  --receipt-output audit/uniqueness-receipt.json \
  --chunk-records 5000000 \
  --auditor-source-commit <40-hex-commit>
```

For the full class, the final binary j-index is 5,729,508,880 bytes including
its sixteen-byte header. It is also the key input to the later external
j-to-node join for the complete degree-5/13/37/59 adjacency pass.

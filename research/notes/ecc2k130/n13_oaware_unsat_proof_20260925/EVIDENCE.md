# Fixed toy UNSAT-proof evidence and replay

Draft [PR #790](https://github.com/aburan28/crypto/pull/790) froze the
[protocol](PROTOCOL.md), source, checker provenance, binary hashes, selected
queries, caps and first pre-outcome checker-control failure at `a10007f`.
Hash-only CI passed before any selected proof generation. The first and only
selected run is [evidence/receipt.json](evidence/receipt.json), from
2026-09-25T15:51:07Z to 15:51:14Z, with `status=success`. No branch or proof
was replaced to select a favorable outcome.

`evidence/proofs/*.drat.gz` holds five deterministic gzip copies of the exact
CaDiCaL text DRAT streams. The receipt records each raw uncompressed byte
count/SHA-256, compressed byte count/SHA-256, the exact #785 query SHA-256,
source and binary hashes, and each cold solver/checker child's command,
wall/CPU, RSS upper bound, UTC interval, exit and stdout/stderr hash. The query
CNFs are regenerated from #781's pinned base plus each fixed assumption; their
bytes agree with #785's archived input hashes. The archive includes the tiny
valid/mutated proof files and checker transcripts, and the failed attempt to
check Q0T0's valid proof against #785's SAT Q0T3 assumption. The first
pre-outcome mutation-test failure is retained separately in
[PREOUTCOME_PREFLIGHT_FAILURE.md](PREOUTCOME_PREFLIGHT_FAILURE.md).

The checker is the vendored `checker/drat-trim.c`, byte-for-byte from the
[upstream DRAT-trim project](https://github.com/marijnheule/drat-trim) at
commit `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`; the upstream MIT
license is retained beside it. The measured Mac checker binary SHA-256 is
`46153180fc81cb4902ef6a648a9de4c775246d63f5c21691d8e5b214a0024326`.
For portable replay, the following command compiles a fresh checker from the
pinned source, verifies the valid/mutated tiny controls, checks every
committed proof against its independently reconstructed CNF, and rejects the
wrong-input control:

```sh
python3 ci_replay.py --evidence evidence
```

The replay does not rerun CaDiCaL or create a new selected proof. The exact
claim boundary and stage-cost ledger are in [RESULT.md](RESULT.md).

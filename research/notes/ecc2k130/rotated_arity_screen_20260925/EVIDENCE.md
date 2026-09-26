# Arity-screen evidence and replay

Draft [PR #778](https://github.com/aburan28/crypto/pull/778) froze the source,
protocol, four arms, resource caps and hashes before the first count. Its
pre-outcome head was `2b565b6`; hash-only CI passed before the first measured
run. The first and only complete run began at 2026-09-25T14:06:47Z and ended
at 14:09:50Z with `status=success` in [the receipt](evidence/receipt.json).
No arm was replaced or rerun to select a favorable count.

The `evidence/raw/m*/producer/` directories contain exact counts and all
8192-row SHA chunk digests; `evidence/raw/m*/verify.json` contains the
independent full replay. The receipt hashes every raw file and each captured
stdout/stderr stream and records all eight command intervals, exits and RSS
upper bounds. The producer used Gray-order x updates, Euclidean inversion and
a precomputed trace mask. The verifier built each x directly from natural
mask bits, used bit-serial multiplication and polynomial-division inversion,
re-derived the trace mask, and checked every row and chunk digest. For 64
SHA-selected masks per arm it additionally reconstructed rational y values,
checked the curve equation and independently doubled twice to the archived
±[4] x key.

Reproduce the archive-only acceptance check from this directory:

```sh
python3 ci_replay.py --evidence evidence
```

This reruns all four independent mask verifiers and compares their
deterministic output to the archived verification JSON. It does not perform a
new producer search or test a PDP solver. The exact per-arm outcomes and their
limited implication are in [RESULT.md](RESULT.md).

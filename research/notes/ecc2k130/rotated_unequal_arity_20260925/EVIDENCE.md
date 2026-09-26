# Unequal-slot census evidence and replay

Draft [PR #784](https://github.com/aburan28/crypto/pull/784) froze the
[protocol](PROTOCOL.md), four arms, caps and source/input hashes at `c798909`
before measurement; hash-only CI passed at that head. The parent balanced
[PR #778](https://github.com/aburan28/crypto/pull/778) had merged, and its
freeze and receipt hashes were pinned in [INPUT.json](INPUT.json). The first
and only measured four-arm run began at 2026-09-25T14:36:44Z and ended at
14:54:05Z with `status=success` in [the receipt](evidence/receipt.json). No
arm was replaced or rerun to select a favorable count.

The `evidence/raw/m*/producer/{low,high}/` directories contain exact counts,
full natural-row SHA-256 digests, signed-column-set digests and every 8,192-row
chunk digest. The arm `producer/result.json` records the rank, tuple product,
threshold outcome, operation counts and normalized nesting check. Each
`verify.json` records a separate all-mask reconstruction, 64 sampled ordinals
per low/high space and source hashes. The receipt hashes every raw and
stdout/stderr file and records all eight cold-child commands, UTC intervals,
exits and RSS upper bounds.

The producer used Gray-order x updates, Euclidean inversion and a precomputed
trace mask. The verifier built every x directly from natural coefficient bits,
used bit-serial multiplication and polynomial-division inversion, rederived
the trace from x, checked each row/chunk digest and the exact column sets,
and checked their containment. For selected rational lifts it independently
constructed y, checked the curve equation and doubled twice to the saved
signed ±[4] x key. Both programs excluded x=1 explicitly; each space has
zero such masks. All four low F/C counts reproduce the independently replayed
balanced parent census.

From this directory, repeat the archive-only acceptance check with a Python
interpreter supporting the frozen source (Python 3.13 was used for the run):

```sh
python3 ci_replay.py --evidence evidence
```

That command checks the frozen hashes and reruns all four independent
low/high-mask verifiers against the archived producer data. It is not a new
producer search, target-support experiment or implicit PDP solver. The exact
counts, gate decision and limitations are in [RESULT.md](RESULT.md).

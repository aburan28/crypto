# Native F5 preparation gate

This is the native replacement for verification of the retained F5 ordinary
preparation, governed by [PROTOCOL.md](PROTOCOL.md). It reads the unchanged
[certificate](../prepared-ic-state-v1/f5-preparation.json) and reconstructs its
matrix independently of the production IC solver. The original 216-query
history and Python controller provenance remain intact.

The [local physical macOS ARM64 replay](local-macos-arm64-replay.json) passes:
`PASS_NATIVE_RETAINED_F5_PREPARATION_REPLAY`. Its source digest is
`3b747c4e409ddd6915aa414e4adb1029beeb8d61282456dba71aff69101cecc5`;
the original certificate file bytes have SHA-256
`15978bbaccf0c7bbf461606166b839eed9abe561a11a62c84156e80a5a3fca44`.
The whole canonical certificate seal differs from the raw file hash by design.
All five F5 corruption/reconstruction controls pass. The complete local harness
also passes 56 controller tests, five worker tests and six integrations; two
legacy integration checks require CI's Valgrind/archive setup. Linux and macOS
CI rebuild the checker and retain their own replay receipts before admission.

| Gate | Retained requirement |
| --- | --- |
| Ordinary query chronology | All 216 original attempts |
| Witness relations | 61 independently checked group relations |
| Negative queries | 155 independently proved geometric negatives |
| Final matrix | Rank 29, 32 dependent rows, zero duplicates |
| Geometry | 63 geometric points, 62 usable images, 29 folded columns |
| Column logarithms | All 29 recovered from the matrix and scalar-verified |
| New queries or target solves | Zero |
| Online timing or speedup | Unknown; no measurement |

The reader preserves `witness` and `proved_unsat` as the original outcomes.
Both complete seals and independent group/matrix checks are mandatory; passing
a byte hash alone is not mathematical verification. No archive extraction,
solver invocation, planted decomposition or target generation occurs.

The [first local verifier failure](first-local-verifier-failure.txt) retains
an input encoding mismatch: archived input point coordinates are decimal
strings, whereas the canonical state uses integers. The corrected reader
decodes coordinates through the checked curve API. It changes no certificate
or historical result.

Use the native command with a new output path:

```sh
/tmp/native-f5-busy busy -- target/debug/icprog f5-preparation-replay \
  --root . --out /tmp/native-f5-preparation-new-replay.json
```

This gate does not establish complete native F5 runtime admission, new natural
yield or fresh-target performance. The complete native controller needs a new
source-bound before-execution registration. All consumed historical controls
stay closed. A later paired comparison uses ecbench after new ordinary-yield,
exposure, reference and calibration gates are satisfied.

# F6 support-local shared-scratch row construction

**Decision: reject for promotion; retain as an opt-in pilot.** This change
acts on the inherited F6-IC support-local Macaulay builder, the actual
matrix-build route on the frozen n17 workloads. The opt-in stream reuses
product and parity scratch across generators and one multiplier schedule
per support-mask/degree-gap pair, while preserving sorted Boolean parity
and exact row order. All 8/8 one-target runs completed and independently
replayed the scalar. Every paired run matched attempts, reductions,
geometric additions, matrix rows/columns, word-XOR counts, and layout
counters. The two T7 build and complete-call improvements fell short of
the preregistered 20%/15% thresholds, and one T1 complete-call run
regressed beyond 5%.

The [protocol](PROTOCOL.md) was committed as `1927314c0` before
implementation or timing. The source and passing unit/build logs were
committed as `b59380310`, and the exact candidate/input freeze as
`f3a060551` before timing. The measured binary SHA-256 is in
[BUILD_IDENTITY.txt](BUILD_IDENTITY.txt). The [IC1 candidate manifests](candidates/),
[input and candidate hashes](FREEZE.tsv), [source hashes](SOURCE_SHA256SUMS),
[all raw run files](runs/), and [raw file hashes](RUN_FILE_SHA256SUMS)
are retained.

The paired workload fixes the prepared n17 Koblitz curve, 62 actual
usable factor-base points, 29 folded columns, three summands, degree
three, archived public T1/T7 targets, imported certified logs, default
algorithm environment, and one Rayon thread. The one-target online
interval starts after reusable setup and includes query, all target PDP
attempts, relation check, descent, and scalar replay. The five exclusive
phase costs sum exactly to each reported online wall interval. No rho
solve was run as part of this stage diagnostic.

| Target | Rep | Original online ms | Stream online ms | Original/stream | Original build ms | Stream build ms | Build ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 3.193 | 3.003 | 1.063× | 2.618 | 2.476 | 1.057× |
| T1 | 2 | 3.475 | 4.290 | 0.810× | 2.894 | 3.428 | 0.844× |
| T7 | 1 | 184.300 | 162.143 | 1.137× | 150.385 | 133.531 | 1.126× |
| T7 | 2 | 157.282 | 150.419 | 1.046× | 129.741 | 124.466 | 1.042× |

The native test compared every support-local row, sorted column and
packed matrix on mixed support masks, colliding Boolean products, zero
rows, and a row-cap failure. Existing direct-packing, fused-packing, and
F6 geometric closure controls passed; the release worker built offline.
The complete target runs additionally fixed all counters and verified
the target scalar. This was an **unisolated macOS host**; its physical
CPU model was not available through `sysctl`, peak RSS is unknown, and
these wall-time ratios are exploratory. They cannot establish a
controlled F6 gain or n83/end-to-end IC gain.

The remaining dominant work is row-product sorting and exact column
collection/packing within this support-local builder. A next candidate
must remove one of those costs while preserving cancellation, exact
column support and the full one-target result. The requested 2×
end-to-end IC improvement remains unestablished, and n83 ordinary
relation yield and same-target rho are still separate gates.

After measurement, Linux generic-admission CI exposed that a new
default-false worker configuration field changed serialized
`effective_config` for older jobs. The PR head now omits the three
new pilot fields when false; explicit opt-in true values remain visible.
This is a post-measurement compatibility fix. All timings above remain
tied to the frozen earlier worker and candidate IDs.

Run `sh research/f6_ic_support_local_stream_20261006/derive.sh` to
rebuild [measurement rows](measurements.jsonl), the [paired derivation
check](DERIVATION_CHECK.json), and [derived hashes](DERIVED_SHA256SUMS)
from the raw runs.

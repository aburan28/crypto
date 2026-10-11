# Exact ranked bitmap columns in the F6 support-local builder

**Decision: reject for promotion.** The opt-in collector assigned an
exact combinatorial rank to each already parity-cancelled Boolean row
term, marked a compact bitmap, and sorted only first-seen column masks
into the original degree reverse lexicographic order. All 8/8 frozen
one-target runs completed and independently replayed the scalar. Every
paired run matched attempts, reductions, geometric additions, matrix
rows/columns, word-XOR counts, and layout counters. On T7, both matrix
build and complete online time were essentially unchanged, far below
the preregistered 20%/15% retention thresholds. The default remains
the original collector.

The [protocol](PROTOCOL.md) was committed as `74dada412` before
implementation or timing. The implementation, exactness tests and
release build were committed as `404966a6b`, and the exact IC1
candidate/input freeze as `f14b07a23` before timing. The measured
[binary identity](BUILD_IDENTITY.txt), [candidate manifests](candidates/),
[input and candidate hashes](FREEZE.tsv), [source hashes](SOURCE_SHA256SUMS),
[raw runs](runs/), and [raw file hashes](RUN_FILE_SHA256SUMS) are
retained. The worker's new default-false flag is omitted from legacy
`effective_config`; an explicit true value remains visible, and its
focused serialization test passed.

The paired workload fixes the prepared n17 Koblitz curve, 62 actual
usable factor-base points, 29 folded columns, archived public T1/T7
targets, imported certified logs, three summands, degree three, default
algorithm environment, and one Rayon thread. The one-target online
interval starts after reusable setup and includes query, all target PDP
attempts, relation check, descent, and scalar replay. Its five
exclusive phase costs sum exactly to every online wall interval. No
paired rho solve was run in this stage diagnostic.

| Target | Rep | Original online ms | Bitmap online ms | Original/bitmap | Original build ms | Bitmap build ms | Build ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 3.021 | 2.931 | 1.031× | 2.418 | 2.360 | 1.025× |
| T1 | 2 | 3.177 | 2.782 | 1.142× | 2.531 | 2.268 | 1.116× |
| T7 | 1 | 152.536 | 152.376 | 1.001× | 127.058 | 126.806 | 1.002× |
| T7 | 2 | 150.203 | 151.674 | 0.990× | 125.113 | 125.431 | 0.997× |

The exactness unit compared columns over all square-free monomials of
degree up to four in small Boolean systems, duplicates, zero rows and
an unsupported-degree fallback. It also compared full packed matrices
on a mixed-support system with product collisions. Existing support-local
streaming, fused-packing and F6 closure controls passed. This removes
flatten-sort/dedup from the opt-in support-local column collector but
did not produce a measurable T7 complete-call gain on this host.

The host was an **unisolated macOS host**; its physical CPU model was
unavailable through the sandbox and peak RSS was not measured. CPU
ratios are exploratory. No n83 ordinary relation, 2× F6 gain, or
single-target IC-versus-rho speedup follows from this stage study.
The next matrix-build probe should time row-product generation and
packing separately before another implementation change. n83 also
requires a higher-arity relation method that constrains partial sums
before enumeration; this n17 three-summand optimization cannot supply
that missing relation yield.

Run `sh research/f6_ic_bitmap_columns_20261006/derive.sh` to rebuild
[measurement rows](measurements.jsonl), the [paired derivation check](DERIVATION_CHECK.json),
and [derived hashes](DERIVED_SHA256SUMS) from the raw run directory.

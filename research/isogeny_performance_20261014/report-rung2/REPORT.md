# BMSS constructor optimization — 2026-10-14

Measured generic Field calls, synthetic 224-bit fixtures, prime degree 101. Counts exclude fixture preparation and verification. Inversions are opaque; this is not an instruction or wall-time ratio.

| Algorithm | Baseline mul calls | Candidate mul calls | Call quotient | Baseline inv calls | Candidate inv calls |
|---|---:|---:|---:|---:|---:|
| linear_algebra | 273659435 | 5640647 | 48.516 | 611 | 209 |
| stark1972 | 10853630 | 10856762 | 1.000 | 1151 | 109 |
| atkin1992 | 547470 | 548379 | 0.998 | 409 | 109 |
| atkin_modcomp | 627414 | 628032 | 0.999 | 417 | 213 |
| elkies1992 | 20109 | 20565 | 0.978 | 255 | 106 |
| elkies1998 | 11127 | 11277 | 0.987 | 249 | 200 |
| fast_elkies | 46958 | 47108 | 0.997 | 108 | 59 |
| fast_elkies_prime | 1073345 | 1073345 | 1.000 | 110 | 110 |

Source: comparison SHA256 `da9c0dd2e21396ad011cdeb4e23bfeb452f150df859689f2e1fb2b9ae157884d`; baseline run SHA256 `8f18f33a0c688125bcdce22ad134dc1f25a38f8b2e9b8d84dbc904ffbf84d64e`; candidate run SHA256 `334faee0d8cec11dc1371eda426a372d6d7396988c9b8142328f7014dacc9ebe`.

The linear-algebra rung moves row-independent series products into a column loop, without retaining an extra polynomial cache. The second rung batches fixed series/Newton denominators above their small-size thresholds. Separate rung receipts distinguish them. A rise in multiplications for batching is reported beside the decrease in inversions.

Native macOS runs used an Apple M4 Pro (14 logical cores, 48 GiB); this shared host had unrelated CPU jobs. No native wall samples were collected or promoted. The four-architecture GitHub matrix is implemented; its receipts must be examined before any other CPU measurement is credited.

Remaining requirements: qualified per-CPU A/A and interleaved A/B timings, additional constructors/finders/path families in INVENTORY.md, named P192/P224 whole-search/replay totals, and actual ReleaseMe asset publication. The earlier 2,000-map search requirement remains open.

Canonical graphs reviewed: docs/index-calculus-scoreboard.html, docs/ic/progress-timeline.json, docs/performance-gains/data.json and comparisons.csv, docs/curves/ and figures/. They describe DLP/rho, unrelated polynomial specialization or curve edges; this constructor operation-count result changes none of those records. It appears in the study graph and benchmarks index.

![Count comparison](COUNTS.svg)


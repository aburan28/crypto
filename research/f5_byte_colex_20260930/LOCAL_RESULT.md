# Byte-colex F5 row packing: local screen

The [protocol](PROTOCOL.md) and baseline build record preceded the opt-in
candidate in commit `02c0f2c98c63764cdb3b46d48c763feb1ebef40b`. The
exhaustive colex-index release test passed; the F5 tests passed 10/10 and
GF(2) elimination tests 6/6. The [two-seed structural receipt](structural_02c0f2c98.json)
contains four complete processes, all seven F5 cases per process, with exact
returned-row fingerprints, canonical row space, term counts, rank, matrix
shape, criterion work and reduction word operations. Its SHA-256 is
`c75630d8ba59bebb78d85568816260cdaf30cf38bcab9e47483d19ec103e81f9`.

The [local paired receipt](local_pairs_adebb63d5.json) contains 44 complete
processes on an Apple M4 Pro: one warmup per arm, five A/A pairs, and five
alternating A/B pairs for each seed. Every output and route check passed.
`tools/isolated_bench.py` cannot reserve a CPU on this macOS host, so these
ratios are **exploratory**. They do not qualify a speedup claim or replace the
required Linux x86-64 isolated pairs.

| Seed XOR | A/A full-call range | A/B full-call median | A/B build median |
| --- | ---: | ---: | ---: |
| `0` | 0.932–1.019 | 1.041 | 1.275 |
| `badc0de1` | 0.974–1.497 | 1.033 | 1.234 |

All ratios are reference/candidate for `f5_n24_m24_d4` at one Rayon thread;
the full-call column uses the benchmark's F5 call interval. The large A/A
spread on the second seed and its 0.695 A/B outlier are retained in the raw
receipt. Neither seed's A/B median is below the frozen 0.95 local stop rule,
so the candidate advances to qualified Linux x86-64 one- and two-thread
screening. The local receipt SHA-256 is
`3d1875573ffa2ad624afc86248fe262d7f1c2ce4fdbc9627c47fb41616c6a7ca`.

This is a Boolean matrix-F5 solver-stage diagnostic, not a one-target IC DLP
or paired rho measurement.

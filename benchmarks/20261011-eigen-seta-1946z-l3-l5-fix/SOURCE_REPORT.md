# L3/L5 round trips after scalar-recovery fix

Protocol frozen before execution on 2026-10-11 19:46 UTC. This run checks
one full KeyGen/Encaps/Decaps round trip each for L3 and L5 of the
secret-degree Eigen-Séta v2 crate imported from `cd47963f9`, with the
local `smul_ct` correction for `[k+1]P = O`, on workspace base
`7f7672d01`. Run `cargo run --release --example run -- L3 1` and
`cargo run --release --example run -- L5 1` in `eigen_seta/`, concurrently.
The example uses `OsRng`; random seeds are not captured.

Acceptance for each arm: exit 0 and print `RESULT 1/1 decapsulations
correct`. Record the reported keygen, encapsulation and decapsulation wall
seconds, and any failure/timeout. These are single correctness samples on
different parameter sets, not a performance comparison, a security
estimate or a constant-time validation.

Results: both exit 0 and print 1/1 correct. On an AMD EPYC 7763 shared cloud
VM with five logical CPUs and rustc 1.99.0 release build, the measured
wall-clock samples were:

| Set | KeyGen (s) | Encaps (s) | Decaps (s) | pk = ct (bytes) | secret f (bits) |
|---|---:|---:|---:|---:|---:|
| L3 | 112.612 | 1.3059 | 47.275 | 1992 | 409 |
| L5 | 259.034 | 1.9862 | 107.221 | 2632 | 536 |

These are one sample per operation per set, run concurrently without
isolation; no spread or uncertainty interval was recorded. The runs preceded
the later public-input validation edit. Exact output and commands are in
`records/run_output.txt`; a quantitative log-time chart is `timings.svg`,
embedded in `report.pdf` from `report.typ`. This is correctness evidence,
not a security estimate, formal constant-time check or comparative benchmark.

# n=83 F6 sign-folded x-only query: retained component

The [preregistered protocol](PROTOCOL.md) compared the retained
sign-folded four-summand pair index with an x-only residual query on the
same K0 n=83 curve, cofactor-projected bases, and public T001 point.
The candidate computes the x-coordinates of both target residuals in
one inversion batch. It calculates a complete curve point only after
an x key matches the pair-sum index, then checks its sign and replays
any four-point witness. The index contents and build path are unchanged.

All six native processes exited 0 and returned the same exact
no-witness result. Every full-base run used 4,054 actual usable points,
8,219,485 unordered pairs, and 4,108,723 signed-sum representatives.
The [status table](status.tsv), [runner](run.sh), six JSONL rows, and
six empty stderr files preserve each run. The final-source
[correctness log](correctness.log.gz) passed the x-key/reference group-law
and exhaustive small-base search controls; the
[geometry log](geometry_tests.log.gz) passed all six module tests. The
logs are losslessly compressed with `gzip -n` and readable with
`gunzip -c`.

| Usable base points | Metric | Sign-folded baseline | X-only candidate | Baseline / candidate |
| ---: | --- | ---: | ---: | ---: |
| 258 | Index build, three-run median | 51.128 ms | 27.219 ms | 1.879× |
| 258 | Exact query, three-run median | 45.601 ms | 8.312 ms | 5.486× |
| 1,048 | Index build, three-run median | 558.781 ms | 763.074 ms | 0.732× |
| 1,048 | Exact query, three-run median | 445.960 ms | 268.333 ms | 1.662× |
| 4,054 | Index build, two-run median | 8.066 s | 7.606 s | 1.060× |
| 4,054 | Exact query, two-run median | 5.066 s | 3.972 s | 1.276× |
| 4,054 | Maximum process RSS | 1,275,953,152 B | 1,276,051,456 B | ~1.000× |

The full runs were interleaved baseline, candidate, candidate, baseline.
Their exact query times were 5.698, 3.626, 4.317, and 4.434 s,
respectively. The full-base candidate median was **21.6% lower** than the
sign-folded baseline, passing the registered 15% exploratory query gate.
Both small-base query medians improved. The full-base build median and
maximum RSS also met their registered limits; the 1,048-point build
median regressed amid pronounced host variation but that was not a
retention gate. The x-only query is retained as an exact alternate method
on `F6SignedPairIndex`; the full-coordinate method remains available for
reference comparison.

The baseline branch head was
`d9f11431a81315fa2612982fd930b405c73e30d5`;
its sign-folded index source SHA-256 was
`9e5ce9ccc3be6a5427a6782ca5fbc65c0f4d6247037776b74969fba9a0af3ec9`.
The candidate index source SHA-256 was
`3f2ef8bd0f2bb07dde0bdb7a9ad7f3da76e0ff1924a9ec04dbd200102ae4125f`.
The small and full candidate probe sources were SHA-256
`77a429eb1d9d6aba39cd259f3ccff402ab1d0c9a91ad1d87bbe2a66caea22bcc`
and `7c9a6338e12972d9ee1b0d10122c02e064bb8106c468c8eb3717bbcc82c338c8`.
The protocol and runner were SHA-256
`c18980dd524067349bf84371b70723992bd34b98bdab509cc71054e8084af8e0`
and `5cfa1070899c73b3289e74013659096ad5e871b01514272c499ec0dbd7333766`.
The baseline full and small probe binaries were SHA-256
`f00655d02e543d324ba709fdb20028585ddfad44fbdf0d486f61bb5e5f6e4adf`
and `f4f96350b1676a7a8140614be0e1a144127a904da6f43163a632f5be50e3b96f`.
The candidate full and small binaries were SHA-256
`850287629d28c5e08d00b62ed42887a91f03f87725fb333df0e8884e8a94f820`
and `3d805bff0ce5307df8cfe1f6d618d0ca65957478f153bb72671d5c747712a2e4`.
Both arms used Rust 1.93.1 release builds on arm64 macOS (Darwin 25.6.0).

This shared host has no isolation receipt, and the raw times vary widely.
The ratios are exploratory four-summand stage diagnostics, not controlled
CPU or complete F6 speedups. This exact four-sum query missed T001 in every
run, consistent with the base's uniform-target coverage ceiling of
`4.662e-12`. No natural relation yield, recovered logarithm, complete
F4/F5/F6 comparison, one-target IC online interval, or paired rho run
was measured here.

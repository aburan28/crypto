# n=83 F6 x bitmap prefilter: rejected

The [preregistered protocol](PROTOCOL.md) tested a 2^25-bit exact bitmap
before each x-only residual hash lookup in the sign-folded four-summand
pair index. Every stored affine x key set its low-25-bit position;
infinity bypassed the bitmap. All possible bitmap hits still used the
full hash lookup, sign check, and curve-group witness replay. The
[rejected candidate patch](rejected_candidate.patch) contains the full
implementation and both frozen probes. The final branch restores the
baseline runtime source and removes the candidate probes from `examples/`.

All six native processes exited 0, returned exact `no_witness`, and used
the same public T001 point. The [status table](status.tsv), [runner](run.sh),
six JSONL outputs, and six empty stderr files preserve the raw run order.
The [correctness log](correctness.log.gz) passed exhaustive small-base
four-sum search and every-stored-key bitmap membership. The
[geometry log](geometry_tests.log.gz) passed all six module tests. Logs
are losslessly compressed with `gzip -n` and readable with `gunzip -c`.

| Actual base points | Metric | X-only baseline | Bitmap candidate | Baseline / candidate |
| ---: | --- | ---: | ---: | ---: |
| 258 | Index build, three-run median | 25.489 ms | 28.688 ms | 0.888× |
| 258 | Exact query, three-run median | 21.442 ms | 9.846 ms | 2.178× |
| 1,048 | Index build, three-run median | 733.014 ms | 349.039 ms | 2.100× |
| 1,048 | Exact query, three-run median | 307.616 ms | 154.636 ms | 1.989× |
| 4,054 | Index build, two-run median | 9.350 s | 8.749 s | 1.069× |
| 4,054 | Exact query, two-run median | 4.407 s | 4.517 s | 0.976× |
| 4,054 | Maximum process RSS | 1,273,856,000 B | 1,279,852,544 B | 0.995× |

The full-base runs were interleaved baseline, candidate, candidate,
baseline. Their query times were 4.042, 4.373, 4.661, and 4.771 s;
build times were 7.404, 7.875, 9.623, and 11.295 s, respectively.
Every full run used 4,054 distinct usable points, 8,219,485 unordered
pairs, 4,108,723 signed-sum representatives, and 2,027 signed columns.
The candidate full-query median was **2.5% higher**, so it missed the
registered requirement of at least 10% lower query time. The monotone
run-time drift and absent host isolation preclude a reliable claim that
the candidate is intrinsically slower. The small-base query medians
improved and the full-base build and RSS limits passed, but the full-query
gate controls retention. The bitmap candidate was reverted.

The baseline branch head was
`e40451d5fafcf5acabc02e98a339aa9081d2711d`, with source SHA-256
`3f2ef8bd0f2bb07dde0bdb7a9ad7f3da76e0ff1924a9ec04dbd200102ae4125f`.
The tested candidate commit was `f7cc15775`, with source SHA-256
`813f19b59b4cb3ffc11295e5cf332b5d0ae575231e9b16a96502aac3dce76cc7`.
The candidate patch SHA-256 is
`dfdebc6823687a41c80b80c6cbe88fe60e165f0456d7e81b49cf522853157763`.
The small and full candidate probe sources were SHA-256
`d53bd0e674e38818b6648dcefd4dbda2407efae32df80c4aaf3ba67bcea7cc15`
and `66a768c2b1740c1de50991aeeaf28bc4758fb185048f4933c99e18cb64ab6d0b`.
The protocol and runner were SHA-256
`218fe2c1c5db94263a6db07538957646f53e9266619f09954074a8e1c7481700`
and `c696966d6849c5012421e6e09dbfd54f80c1902d2416509c8ea3e2b8014c7089`.
The baseline small and full binaries were SHA-256
`3d805bff0ce5307df8cfe1f6d618d0ca65957478f153bb72671d5c747712a2e4`
and `850287629d28c5e08d00b62ed42887a91f03f87725fb333df0e8884e8a94f820`.
The candidate small and full binaries were SHA-256
`3e67794451a970b724576a864e4cdbd39c0aa758d52ba9305b384c87daae3d7d`
and `4d48bea301c8c9d12dc8f23c04adf16c0dc1bde74c6bc1835e3e48e8647f4b90`.
The [build log](build.log.gz) records the offline release build with
Rust 1.93.1 on arm64 macOS. Compiler temporaries were redirected to
the SSD after the system volume fell below 400 MiB free; the initial
interrupted test build produced no timing row.

These are exploratory four-summand stage timings on a contended host,
without a CPU-isolation receipt. No controlled speedup follows from
the small-base numbers. This K0 base's uniform-target four-summand
coverage ceiling is `4.662e-12`, and T001 was an exact miss in every
run. No ordinary relation yield, recovered logarithm, complete F6/F4/F5
comparison, one-target IC online interval, or paired rho solve was
measured here.

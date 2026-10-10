# C2800 early root-slot prefetch

This variant issues the root-table home-slot prefetch immediately after each
S3 partner key is canonicalized. It gives the memory fetch the remaining root
arithmetic in the eight-state batch to overlap, while the later lookup pass
checks witnesses in the original order. The frozen N53 T025 control recovered
the expected scalar, matched every non-timing field of the prior late-prefetch
report, reached rank 2800, and used 24.02 GiB peak RSS.

The source is
`examples/s3_partial_prepared_fullrank_c2800earlyprefetch_targetspan.rs`
at commit `188f4007b`, SHA-256
`ac8453008ef252c00651499a81a6c6e2942977896d1cc425c7f007786ace12f5`.
Its exact candidate ID is
`IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h0885f0c36bda`.
The control freeze binds its binary, T025 public point, previous report,
runner, and 32 GiB memory envelope. The raw receipt and compressed report are
under `local/verified_t025/`.

| T025 same-point diagnostic | Late prefetch | Early prefetch |
| --- | ---: | ---: |
| Verified scalar | 429728658711 | 429728658711 |
| First witness point indices | `[55035, 67241, 57152, 189185]` | same |
| S3 calls | 512 | 512 |
| Online ms | 0.310 | 1.828 |
| Charged target PDP ms | 0.161 | 1.668 |
| Peak RSS GiB | 24.08 | 24.02 |

The two intervals are separate local executions. The same late-prefetch
binary recorded 1.818 and 0.310 ms online in its two T025 controls, and the
eight-point late-prefetch screen failed its reusable-setup noise check.
Therefore these intervals are stage diagnostics, not a controlled ranking of
the two prefetch schedules. The next experiment should compare lookup modes
inside one prepared index with order balance and independently verified
relations, then carry a promising mode to an isolated one-target panel.

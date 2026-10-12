# PS0 PDP oracle baseline, five rungs, A/A and counted operations

[Source study](../../research/pdp_oracle_speedup_20261011/baseline_v1/README.md) and [run report](../../research/pdp_oracle_speedup_20261011/reports/run-20261011T2107Z/report.pdf) remain the source of truth. These browsing copies were added during integration on 2026-10-12. The PDFs and SVGs are byte-identical copies; each PNG renders its copied SVG at 150 dpi.

Baseline medians 337.7 / 375.4 / 393.3 / 421.4 / 452.3 ns per pair; largest A/A spread 8.81%; accounting result

Host and scope: 4-vCPU Intel Xeon 2.80 GHz, shared Linux container; CPU 3 pinned; 0/50 timed runs marked contended; no certified isolation level. The integration adds no measurement or independent replay.

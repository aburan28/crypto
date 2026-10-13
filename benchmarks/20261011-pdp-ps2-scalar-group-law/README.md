# Scalar group-law PDP kernel versus baseline, five rungs

[Source study](../../research/pdp_oracle_speedup_20261011/ps2_scalar/README.md) and [run report](../../research/pdp_oracle_speedup_20261011/reports/run-20261011T2255Z/report.pdf) remain the source of truth. These browsing copies were added during integration on 2026-10-12. The PDFs and SVGs are byte-identical copies; each PNG renders its copied SVG at 150 dpi.

29.8–40.8 ns per own pair; group-scalar/baseline per-pair ratios 0.0668–0.1225; 400 corpus root sets agree; library promotion and Gray update remain partial

Host and scope: 4-vCPU Intel Xeon 2.80 GHz, shared Linux container; CPU 3 pinned; 6/125 timed runs marked contended, retained with sensitivity table; no certified isolation level. The integration adds no measurement or independent replay.

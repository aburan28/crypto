# Prime rootfinding collector buffer reuse

The report and chart here are byte-for-byte copies of the visual results in [`research/ecbench_rootfinding_alloc_20261014`](../../research/ecbench_rootfinding_alloc_20261014/). The study directory contains the frozen spec and plan, 384 sealed executions, 288 measured replays, exact cross-session comparisons, raw failure and correction records, and the rendering source.

![Paired L0 wall-time diagnostics](timing-ratios.svg)

The candidate's batch-8 relation-phase wall ratio to the first baseline was 0.937 (95% target-cluster interval 0.786–1.122); whole-solve wall ratio was 0.993 (0.803–1.251). Both intervals include parity, and the repeated baseline plus rho A/A controls show substantial shared-host variation. No qualified speed claim follows. See the [one-page PDF](REPORT.pdf) and [full study report](../../research/ecbench_rootfinding_alloc_20261014/REPORT.md).

SHA-256: `REPORT.pdf` `a5424c5750fa46daefa4b60f5c88067fc93200b227e8bdb7aea3c0601c6d2645`; `timing-ratios.svg` `d3b4e6e6db30eac77bede84cb12716f0fd0bf57950608b9587dbf93a8393d7a5`; `timing-ratios.png` `57f9091f4c235d924c50c36553c7dcf9c4f17c54934bd22bc3f3f01027a9316d`.

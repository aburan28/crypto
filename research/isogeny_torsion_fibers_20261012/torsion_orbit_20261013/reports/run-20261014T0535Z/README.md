# Run report: exact torsion-orbit support screen

UTC run start: 2026-10-14 05:35. Protocol commit `1d01a39a25`.
Host: macOS Darwin arm64, shared/L0. The enumeration is deterministic;
no random seed, wall-time comparison, or confidence interval applies.
Source, binary, input and raw hashes, commands and host are in
[RUN.md](../../RUN.md). The native program requested and completed twelve
rows: three factor-base arms at each of four subgroup orders. There were
no runtime failures. The first compile attempt failed before execution and
is retained as [compile-v1.log](../../compile-v1.log).

The 36-lift arm and 12-point arm had exactly equal two-sum support:
64/103, 74/439, 75/1613 and 78/6427. The 36-distinct arm reached
103/103, 351/439, 559/1613 and 642/6427. These are exact finite-group
fractions, with no sampling interval. The [visual](../../graph-support.svg)
and [run PDF](report.pdf) display those counts. The [result report](../../REPORT.md)
explains the proof, interpretation and next gate.

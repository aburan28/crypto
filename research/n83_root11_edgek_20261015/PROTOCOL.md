# Root11 wide3 factor-base capacity edge at N83

This fresh-target screen tests K3500 and K3550 against the independently
replayed K3450/root11/wide3 baseline. The three pipelines use the same
binary source SHA-256
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`,
symmetric S3, sparse canonical root selector, fixed five-byte root index,
wide3 Bloom prefilter, rank seed 29, four rank threads, and one target
thread. Every process starts with an empty base, root index, and factor logs.
K is the folded column count; the actual usable factor-base size is
`166*K` only after independent base replay.

`generate.py` fixes 1,000 public points before the six-arm allocation from
domain `n83-root11-edgek-fresh-1000-20261015-v1`. The public-point SHA-256
is `c75cbb92ef1a23057272cf3481cc649b12495135c75affd1c2610205a42b7aef`;
the scalar sidecar SHA-256 is
`ca95d57ff091b8fc7de8aa2fd22120cadac94a7c7a697e3c757aebc53115a159`.
The generator excludes all earlier N83 panels, the 4,000-point final
tournament, and the five single-target Tier 1 points. The scalar sidecar
stays outside every producer input and timing interval.

Use one Ubuntu 24.04 cloud allocation exposing at least four CPUs and at
least 14 GiB observed peak RSS capacity plus 256 MiB host headroom. Require
2 GiB free disk. Run K3450, K3500, K3550, K3550, K3500, K3450 in separate
fresh processes, each with a 1,200-second child cap and 14 GiB observed
RSS guard. Preserve CPU affinity, cgroup CPU/memory limits, throttle and
OOM deltas, host topology, source/binary hashes, stderr, exact base/index
counts, target phases, setup costs, and every timeout or failure. This is
a cloud-instance exploratory secondary many-point factor-base screen.

Independently replay each factor base, target scalar, and four-point
relation. For each target, average the two verified online intervals of a
pipeline. Use 40 contiguous blocks of 25 target indices and 10,000 paired
bootstrap resamples at seed `2026101580+i*4+j` for ordered candidate
positions `i,j`. Select the lowest arithmetic mean only if the upper
endpoint of its 99% interval against both others is below zero. Otherwise
record `BOUNDED_UNKNOWN`; any incomplete or ineligible arm also leaves the
choice unknown. Preserve all six arm rows and the complete comparison.
The selected or descriptively lowest eligible pipeline can enter a later
fresh-target comparison under its own frozen admission rule; this screen
does not alter the four-candidate tournament already running.

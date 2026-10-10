# N53 C1200 one-versus-four target-worker screen

This experiment keeps the same exact N53 C1200 factor base, root index,
relation collector, and complete rank-1200 precomputation. The frozen source
and release binary are identical in both arms. The reference uses one target
query worker; the parallel candidate starts four independent rotating target
searches and stops them cooperatively at the first verified relation. The
prepared index and base-log system remain target-independent. Each solver
process receives one public point, and all target workers' completed work is
charged inside that point's online interval.

`SCREEN_FREEZE.json` binds a correctness calibration pair on T025, followed
by eight previously unused public points T181–T188. The T025 pair is excluded
from the held-out ratio. The run order alternates AB/BA, with one repetition
per arm and point, 14 preparation workers, at most four target workers, a
32 GiB RSS envelope, and a 120-second process timeout. A fresh helper process
provides per-solve RSS and resource counters. The runner preserves failures
and stops at the first invalid or over-budget row. After the T025 pair, it
pauses at a frozen verified prefix so the report and resource gates can be
checked before the held-out panel starts.

The predeclared local advancement threshold is an eight-target geometric
mean of at least **2×** for one-worker/four-worker online wall time, with all
16 held-out scalars and phase sums verified. Target S3 calls, state probes,
preparation times, and phase costs remain diagnostics. The same physical CPU
resources are available to both arms; this ordinary macOS host has no
exclusive CPU partition or thread-affinity receipt. Any wall-time result is
therefore exploratory under the repository's CPU isolation rule.

# Host memory limit implementation

Added before execution. This macOS host rejects both `ulimit -v 8388608`
and `ulimit -d 8388608` with `setrlimit failed: invalid argument`.
The runner therefore samples the native process's RSS every 50 ms and
terminates it at 7 GiB sampled RSS, leaving a 1 GiB guard below the
registered 8 GiB envelope. It also stops after 120 seconds. This is a
monitored limit, not a kernel-enforced hard limit: an allocation spike
between samples could exceed 8 GiB. Preserve the maximum sampled RSS and
`getrusage` peak RSS, and classify any such overshoot as an envelope
failure. The primary outcome is algebraic feasibility and exactness,
never an isolated wall-time speedup.

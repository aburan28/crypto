# Run-wrapper amendment before any 75% arm result

The frozen first attempt began at 2026-10-13T11:14:42Z. Its K1798 producer
wrote a complete one-target summary and target row, but macOS
`/usr/bin/time -l` returned status 1 after the producer because sandboxed
`sysctl kern.clockrate` was denied. The initial `runs/primary/k1798_dense`
row retains its exact output and exit 1; it is an infrastructure failure,
not a promoted measurement. The following K2000 and K2202 attempts were
interrupted while stopping that script. Their partial files and exit 1 rows
remain in `runs/primary/`.

Before obtaining any K2000/K2202 results, `run.sh` changed only the outer
timing wrapper to `/usr/bin/time -p`, which has a successful local zero-exit
preflight. The producer source/binary, public points, factor-base arms,
order, 900-second cap, rank seed, thread counts, and target timing boundary
are unchanged. The corrected primary rows go to write-once
`runs/primary_v2/`. The producer summary contains the target's exclusive
online phases and allocated table bytes; `time -p` does not report peak RSS.
This primary attempt therefore leaves its peak-memory field missing.

Before starting the secondary panel, `run_panel_v3.sh` pins a checked
direct-child wrapper (`run_with_usage.py`) and writes one kernel-reported
peak-RSS receipt per arm in `runs/panel/`. It retains the same 900-second cap,
source/binary/input hashes and reverse arm order. The wrapper's preflight
receipt is `usage_preflight.json`. The independent checker accepts only
complete zero-exit rows and checks each panel resource receipt.
The panel starts its producer with a clean environment containing the declared
KIC flags, so optional Bloom and sparse-canonical paths are explicitly off.

## K2202 local resource guard

The K2202 primary attempt produced no target or summary row. After about
9.5 minutes, the producer was in macOS uninterruptible I/O state with a few
MiB resident and almost no CPU progress; `vm.swapusage` reported all 13,056
MiB used. The raw process and swap snapshots are
`k2202_pressure_ps.txt` and `k2202_pressure_swap.txt`. The operator
interrupted this owned process to release memory pressure; the outer script
exited 130 before writing `exit_status.txt`. The write-once partial directory
remains in `runs/primary_v2/k2202_load75`. This is a censored resource row.

`run_panel_v3.sh` was frozen but never launched. Before any panel arm,
`run_panel_v4.sh` removed K2202 under this resource guard and retained
K2000 followed by K1798, reversing their primary order. Both use the same
public 64-point panel, 8 GiB packed table bucket, 900-second cap, source,
binary and clean environment. The panel selection is therefore among these
two completed local candidates only. K2202 requires a host with enough
memory headroom for a separate rerun.

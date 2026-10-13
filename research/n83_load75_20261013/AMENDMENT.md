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
`runs/primary_v2/`; the secondary panel remains `runs/panel/`. The producer
summary contains peak RSS and the target's exclusive online phases. Its
independent checker accepts only complete zero-exit rows.

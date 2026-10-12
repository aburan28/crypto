# Frozen four-thread F5 control after the one-thread 2× gate

The one-thread experiment defined in [PROTOCOL.md](PROTOCOL.md) passed its
2× gate on frozen primary `f5_n24_m24_d4`. This control measures the same
three source modes with four Rayon workers; it does not alter the candidate
or tune based on the one-thread timings. The measured candidate code is
commit `5265f0e5f3be835a6af8a1efd64dfcacf3417ef8`. Only research
documentation and CI workflow wiring have changed since that code head.

Use the same [paired_f5.py](paired_f5.py) driver and flags as the one-thread
run, with `--pairs 5 --threads 4`. Freeze all seven F5 cases and the same
four seed XORs: `0`, `badc0de1`, `5eed2026`, `f5c02a28`. The primary
case remains `f5_n24_m24_d4`. Require a Linux x86-64 host with AVX2,
BMI2 and four allowed CPUs; otherwise record `unsupported_host` or the
driver's CPU-count error before timed calls. Pin the process to the first
four allowed CPUs and set `RAYON_NUM_THREADS=4`. On each seed, warm each
arm once, take five default/default A/A pairs, then five rotating
default/prior/new triads. Keep every process output, status, source and
binary digest, CPU features, affinity, load, phase times, raw and canonical
row fingerprints, row and column counts, rank, pruning and word operations.

`wall_ms` remains the full matrix-F5 call inside the example; input
fixture construction and post-call fingerprinting are outside it. For each
seed and case, report the medians and exact five-pair bootstrap 95%
intervals of default/prior, default/new and prior/new complete-call ratios,
with default A/A ranges. Require the same exact output and smaller-case
checks as the one-thread protocol. A four-thread 2× claim requires the
frozen default/new median **and** lower interval bound above 2.00, all
holdout primary medians above their A/A maxima, and no smaller-case
default/new median below its A/A minimum. A useful four-thread gain short
of 2× requires frozen default/new median above 1.05 and lower interval
bound above 1.00 with those same guards. If neither gate passes, report
the regression or uncertainty and keep this combination opt-in.

Do not rerun or select an eligible host based on timing. The prior
one-thread result stands on its own measured host; if this control lands
on a different CPU, do not divide or compare its absolute milliseconds
with the prior run as a scaling measurement. No default behavior changes
in this control. This is a solver-stage diagnostic, not an IC online-time
or DLP speedup measurement.

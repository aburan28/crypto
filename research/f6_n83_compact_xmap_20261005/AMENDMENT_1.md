# Small-base regression adjudication

Registered after the first panel, before the replay below. The original
full-base gates passed, but its single small-base process reported a
dimension-8 query median of 1.640 ms for baseline and 15.913 ms for the
candidate, despite dimension-10 improving. That 9.70× apparent
dimension-8 regression is too large to ignore on a contended host.
Keep both frozen binaries and all previous observations unchanged.

Run only the existing small-base probe, in baseline, candidate,
candidate, baseline order, with no builds between arms. Each process
contains three dimension-8 and three dimension-10 repetitions on the
same public T001 and has a 120-second cap. Preserve all raw outputs,
stderr and statuses. Compare each process's median dimension-8 query
time, then take the median of the two process medians per arm. If the
candidate is more than 20% slower than baseline on dimension 8,
**reject** the runtime change despite its first full-base gate; the
smaller-case regression risk outranks the full-base benefit. If it
does not exceed 20%, retain only if the original full-base gates and
all correctness tests still pass. Dimension 10 remains a diagnostic.

This replay is an exploratory regression adjudication, not an
isolated-host CPU speedup. Preserve the first panel and failed SSD
test-build log; neither is overwritten.

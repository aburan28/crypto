# Feasibility probe (parallel, superseded)

The first pass of `run.sh` with `JOBS=3` (three runs sharing four cores),
kept as evidence and not used in `RESULTS.md`: the two algebraic arms are
priced by *measured* wall (the runner declines to price by their partial
word-XOR count), so concurrent runs share that price.  The container
restarted before `k0_31 sym4` and `k1_17 x-m3` finished; their logs
are here.  The serial pass in the parent directory is the evidence.

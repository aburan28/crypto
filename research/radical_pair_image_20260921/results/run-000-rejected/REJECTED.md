# Rejected run 000

This run is incomplete and is excluded from every result table.

The runner used `msolve -g 0 -v 2` in an attempt to preserve the F4 protocol
without printing a large Gröbner basis.  In msolve 0.9.5, `-g 0` continues into
solution parametrization instead of stopping after the Gröbner stage.  The
first 27 tiny-cell processes completed, but the next process entered this
different workload.  The run was interrupted rather than mixing it with the
intended solver-stage comparison.

Run 001 changes only the output mode to `-g 1`, which stops after the same
Gröbner computation and writes the compact leading ideal.  Inputs, target
selection, variants, repetitions, timeout, thread count and random seed remain
as frozen in `contract.json`.

No measurement from this directory is accepted or quoted.

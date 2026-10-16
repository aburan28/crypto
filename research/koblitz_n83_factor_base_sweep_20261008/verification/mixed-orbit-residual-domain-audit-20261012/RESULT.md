# Independent residual-domain replay for the mixed-orbit gate

The outer event replay previously checked a candidate's exact binary-curve
point equation and its small-base indices, then accepted the worker's reported
residual points without independently checking their membership in the pinned
K=2,048 envelope's suffix. A group-valid equation can use arbitrary subgroup
points, so the group check alone does not prove that a partial belongs to the
declared two-residual search domain. The worker's own lift restricts points to
that domain, but this is an independent evidence check of the worker output.

The replay now takes both byte-pinned compressed objects. For a candidate it
checks the exact small and envelope point counts, loads the referenced small
rows, scans the envelope from the K=1,182 prefix boundary, and requires every
reported residual point to occur in that suffix before accepting the group
equation or graph row. The outer run already checks the objects' SHA-256 values
before this replay. A candidate-free search keeps its earlier early return;
the suffix scan is charged only when a partial is reported.

The source test includes a positive group-valid partial with residual points
in the suffix and a negative **group-valid** partial with residual points
outside it. The negative case is rejected specifically for suffix membership,
so it cannot be rescued by its correct group equation. Existing coefficient,
counter and completed-row replay checks still pass. The full study Python
suite passed **46/46**; its raw output is [study-python.log](study-python.log).

This is a source-only correctness repair. It launched no retained-N83 search
and adds no experimental charge. Partial yield, graph closure, rank and
complete cold runtime remain unmeasured. Any later search needs the separate
authorization in the [current protocol](../mixed-orbit-search-12600-source-20261012/PROTOCOL.md).

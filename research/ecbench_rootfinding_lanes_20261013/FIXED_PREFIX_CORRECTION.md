# Fixed-prefix constructor correction

The first execution of the frozen [fixed-prefix protocol](FIXED_PREFIX_PROTOCOL.md) stopped at the curve-identity assertion before generating any probe. The [raw first-run log](fixed-prefix-initial-failure.log.gz) records `left: 1046999`, `right: 562427` for the subgroup order and a failed test exit. Its result is an implementation failure, not coverage evidence.

Cause: the protocol and first test used `prime_instance_for(20, 59297)`, which selects the 20-bit bench roster before considering the seed. The sealed ecbench `prime_search` workload uses `find_prime_order_curve(20, 59297)` directly. The latter produces the frozen public curve with `r = 562427` and the target coordinates in the test. The assertion prevented accidental measurement on the wrong curve.

Correction: the test now calls `find_prime_order_curve(20, 59297)`. This replaces only the constructor clause of the frozen protocol. All eight public targets and algorithm seeds, B48 bases, batch-8 walk, 20,000 attempted-point cap, 8,192-surviving-probe prefix, milestones, rank and cost definitions, and failure rule remain as declared. The correction is committed before the next execution. The first log remains unchanged; the rerun will have its own output and source commit.

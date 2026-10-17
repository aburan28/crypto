# Next bounded mechanism: exact nondecreasing source constraints

Apply ONLY to the proved symmetric norm predicate with identical source domains. For each adjacent pair of d-bit sources, process most significant bit first. Start equal=1, less=0. For each bit x,y, set term=equal AND (NOT x) AND y; update less=less XOR term, equal=equal AND NOT(x XOR y). The updates use the old equal. At every step less and equal are disjoint, so the final less XOR equal is exactly x<=y, including equality. Enforce it for all four adjacent pairs.

For fresh independent inputs, an upper bound is3d-2 AND gates per comparator (first-prefix constants simplify), hence208 added AND gates at d18. This is a derived circuit bound, not a native count or speed measurement. Less-first comparisons or componentwise bit ordering are invalid: binary1<=binary2 despite their low bits ordering oppositely.

New pre-registration oracle must exhaust all adjacent pairs for d1..5 and compare the compiled truth table to integer order, then compare complete native SORTED norm-primary models against independently sorted literal roots. Restore every unique source permutation, including zero/repeated sources, and require exact equality to the full original norm-primary set. Preserve r and kernel bit. Original S3 auxiliary restoration remains separate and is never inferred.

Proposed initial hard prediction, conditional on the current limit2 outcome: with exact ordering constraints, find a literal-replayed norm witness or complete native norm exhaustion in60seconds under4GiB/CPU0 and90-second outer cap. Freeze a new source/plan, retain UNKNOWN as unknown and never divide by120 uniformly. This is a changed-search hypothesis, not an increased-time retry or qualified all-roots speed claim.

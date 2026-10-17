# Cross-session checker correction

The first `compare_sessions.py` invocation saved `candidate-vs-baseline-a0.json` and `baseline-drift.json` with apparent `phases` mismatches for the IC arms and exited 1. Inspection of the first batch-8 record showed that the differing phase keys were exclusively `wall_ns` (for example, baseline/candidate relation phase 1,050,207/517,922 ns on the first warmup record). Phase operation counts and native counters were identical in that inspected record.

The checker was corrected to remove only `wall_ns` from each phase for deterministic equality. It keeps the original failed outputs and writes corrected comparisons to new `*-corrected.json` files. The separate timing analysis still uses wall times. This is a checker classification error, not evidence of a changed relation stream or a measured improvement.

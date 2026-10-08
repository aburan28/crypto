# Follow-up: larger cache with the selected register schedule

Status: preregistered before compilation or timing.

The first screen established that B32 cache8 improves the occupancy-matched
cache4/min2 control by 1.75% (95% CI 1.0085..1.0266), but the min2 compiler
schedule remains slower than selected.  This follow-up compiles B24 cache6 and
B32 cache8 with `ECC_MINBLOCKS=3`.  The 102,400-byte SM shared-memory capacity
still limits both candidates to two resident blocks; the launch bound asks
ptxas to retain the selected 80-register schedule instead of the 126/128
register min2 schedule.

The matched controls are selected B16/cache4/min3 and the existing
B24/cache4/min3 and B32/cache4/min3 binaries.  No other flag changes.

Before timing, both candidates must compile, pass byte-identical built-ins,
normalized state and DP comparison at run IDs 0 and 139, the actual-walk
shared-X unit, and memcheck/initcheck/synccheck.  Timing uses three interleaved
equal-work repetitions of exactly 34,359,607,296 scalar updates per arm.  A
candidate qualifies only if the paired 95% confidence interval over selected
is wholly above 1; promotion still requires fresh held-out confirmation.

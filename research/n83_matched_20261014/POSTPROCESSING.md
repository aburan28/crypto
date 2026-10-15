# Post-campaign verification revisions

The native campaign completed all32 declared processes with exit0. The CPU3
reservation ended before any independent replay, compilation of post-campaign
tools, analysis or rendering. The original campaign, measured sources,
executables, fixtures and prospective protocol remain unchanged.

The exact frozen v2 checker was run first, with its own native resource meter.
It rejected copied floating-point metadata: its receipt reports227 invalid
metadata/row findings and11 scalar replays. These are checker-attempt failures,
not227 failed solver answers. The failed receipt, stdout, stderr, exit status
and resource record remain in `records/postprocess-20261015T0018Z/`.

The independent native diagnosis found127 discrepancies between original JSON
and copies parsed and serialized by the Rust controller. Every discrepancy is
exactly one decoded f64 ULP; none is an integer, scalar, coordinate, string,
type or key discrepancy. A separately versioned [v3 checker](verification/revisions/20261015-v3/README.md)
allows at most one adjacent finite decoded f64 value only at explicitly named
metadata-copy paths. Original file digests, integers/u128, scalar and witness
arithmetic, seeds, counters, matrix, resource and A/A gates remain exact. The
existing arithmetic-derived eight-epsilon checks are unchanged. Native tests
include two-ULP and one-nanosecond rejection,30 negative cases and442 actual
copy bindings; source and executable identities are recorded.

The separately metered v3 replay passed32 processes,218 canonical scalar
answers (104 IC,114 rho),104 IC four-point witnesses and five identical A/A
work pairs, with zero invalid or missing answers. Its audit cost was
1.1674539 seconds wall,1.089206 seconds user CPU,0.075961 seconds system CPU
and109332 KiB peak RSS, with no observed contention. These costs are outside
both solver arms. This verifies answers and exported witnesses/file bindings;
it does not export or independently validate the complete IC rank/log transcript.

The exact frozen `runner collect` was attempted next. It returned exit2:
`resume raw records do not reproduce the frozen observation`. Its strict Rust
Value comparison encounters the same additional serialization round trip.
No analysis input was created by that failed command. The frozen controller
has not been edited or relabeled as passing.

A separately versioned [collect-only postprocessor](postprocessing/collector-only/README.md)
passed the complete v3 gate and sealed all336 snapshotted files by SHA256 before
extracting32 process rows,192 cohort rows and26 auxiliary rows in the original
analysis-input schema. Its arbitrary-precision JSON copying retains decimal
tokens without another floating-point serialization step. It has no solver
launch operation and changes no observation. Its six native tests passed.
The [extraction receipt](records/postprocess-20261015T0018Z/collector_extraction_receipt.json)
records its source, executable and configuration identities and binds the
failed original collector attempt.

The original frozen native statistical analysis then ran successfully with
`all_required_rows_valid:true` and `overall_primary_winner_eligible:true`.
Its executable remains SHA256
`b1820abad3d17cb4f40f7d4d1a0a6c92ade2dfa056a43a0c4eee2d315c3a8346`.
All declared estimators, bootstrap seeds and experimental logic are unchanged.

This is an explicit post-campaign validation/presentation correction. It is
not a rerun, a change of target panel or method, a different time boundary,
new parameter selection, or an alteration of measured raw output. A future
producer/serializer correction would need a new protocol and new measurements.

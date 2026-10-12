# Frozen replay repair after the one-shot run

This is a narrowly scoped verifier repair, declared before replaying the
archived target data. The original protocol, producer, configuration, target
streams and 16 raw cells are immutable. The only panel run is
[GitHub Actions 36722040881](https://github.com/aburan28/crypto/actions/runs/36722040881)
on merge commit `5548b284c1c5a95538f05e23a2a7a19322300ebc`. Its producer
reported `PASS_PANEL`, but its verifier raised `AttributeError` while
constructing the first unmetered Vélu map, before any raw cell was replayed.
The original `replay.json` is a failed verifier receipt, not a successful
independent check. The original archived panel remains unaccepted until the
repaired verifier finishes.

The cause is the pilot module's process-global replacement of
`oriented_velu.Koblitz` with a metered curve class. The verifier used a bare
field, whose lack of a `meter` attribute made that constructor fail. The
repair restores the bare curve constructor in the verifier process before
geometry is constructed. It adds a direct pre-panel map/dual regression test,
pins the changed verifier/test/workflow source in `FROZEN_REPLAY.json`, and
disables another panel dispatch on the updated default-branch workflow.
The old `FROZEN.json` remains byte-for-byte unchanged and still pins the
original producer and run. The repaired verifier must check the old lock,
the exact result and failed-replay SHA-256 digests, all 16 cell hashes, all
8,192 target cases, map controls, independent rank/scalar recovery, point
law samples, and cost-ledger reconstruction. Its output goes to a new file;
the failed receipt is never overwritten.

This repair does not authorize another target generation, altered base,
changed success criterion or result-dependent tuning. A second verifier
failure remains a `FAIL` receipt and must be investigated openly. A replay
`PASS` accepts only this degree-21 toy panel; it leaves n131 PDP yield,
full-ECDLP cost and rho crossover unset. Preserve both receipts and the
raw archive in the outcome PR and update the scoreboard from that archive.

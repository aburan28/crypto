# Procedure deviations and availability

- The initial full-library test build could not read fourteen fixtures excluded
  by sparse checkout. The source fixtures were materialized unchanged, and the
  build was retried. Both logs are retained in the final evidence.
- Compilation completed between process inspection and the pause command.
  That command named the exited compiler rather than its newly started test
  process. The test process consequently overlapped the initial 29-bit session;
  it was suspended once detected and remained suspended for the subsequent
  session and audits. This violates the protocol's no-test-overlap resource
  condition. The records retain their actual contention/isolation observations.
  No wall-time result is admitted or interpreted. Deterministic counted units
  remain the prespecified primary diagnostic, conditional on exact replay;
  the overlap is not silently removed from the evidence.
- The Conductor CLI/workflow and Sage were unavailable. Conductor check was
  attempted; prior Sage artifacts are historical evidence, not a new Sage run.
- Exact representation manifests were added to the native analyzer after its
  initial source hash was committed, before measurement. Its statistical
  procedure, seed, sample count, inputs and decisions did not change. Source
  and executable hashes for both versions remain separate records.
- ecbench's documented `--allow-busy` mode permits below-L2 operation-count
  diagnostics. Native CPU reservation was attempted with `--cpus auto`; actual
  isolation levels are recorded. No isolated wall-cost advantage is claimed.

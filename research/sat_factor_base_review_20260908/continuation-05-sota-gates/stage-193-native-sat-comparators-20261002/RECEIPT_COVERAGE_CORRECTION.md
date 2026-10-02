# Stage 193 receipt-coverage correction

The solver, WDSat provenance, and CryptoMiniSat provenance receipts were fully
replayed from the first composition. The generic additive charge originally
read every `metrics.json`, but setup, input-copy, Rust build/test, composition,
and verifier receipts outside those named roles were not all replayed.

Before publication, the finalizer was tightened to collect every
`receipt.json`, validate its schema and process accounting, rehash its exact
stdout, stderr, and metrics artifacts, compare its embedded process record to
the metrics file, and require the set of receipt-covered metrics to equal the
set charged by the stage. The previously valid `21/21` result with SHA-256
`9f75dc42d7e1ec5eb2fb25b02266cebd4ae0c3c03bfb9e72869cbf50b5f6d6c8`
is preserved under `development/superseded-receipt-coverage/`. Solver
measurements and interpretations do not change; the preservation, corrected
build/test, composition, and replay are charged additively.


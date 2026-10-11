# P-192 weighted-CM independent verifier

This standalone Rust workspace verifies the frozen Role-1 P-192 weighted-CM
preflight artifacts. It has no dependency on the repository root crate or the
producer package, and it performs no scientific BOX census.

The build requires `P192_WCM_VERIFIER_COMMIT` to be the exact nonzero,
lowercase 40-hex Git `HEAD`. A release build additionally requires a clean
worktree. The runtime preflight command refuses anything other than such a
clean release build:

```text
p192_weighted_cm_verify preflight \
  --source /absolute/path/to/RUN_DIR \
  --out /absolute/path/to/RUN_DIR/independent-verification.json
```

The verifier accepts only the frozen canonical JSON and binary schemas. It
independently rederives the P-192 identity and CM order, replays primality
evidence, regenerates the all-and-only factor base, replays Role-1 controls,
and regenerates all 1,025 REF-0 candidates, dispositions, certificates, and
roots. It writes only `independent-verification.json`; agreement and receipt
artifacts remain supervisor-owned.

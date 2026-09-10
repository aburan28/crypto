# Licensed-host Magma reproduction

This procedure consumes the twenty immutable Stage 13 inputs and runs only the missing Magma lane. It does not regenerate the source systems or reuse planted points as witnesses. The runner first verifies the Stage 13 archive and merge-commit custody, then hashes every manifest, ANF, CNF-XOR, Magma input and source-binding receipt.

## Required tools

- A licensed Magma executable with the direct Faugère F4 implementation.
- A `minisat` executable in the path used by Magma. Magma's documented `SAT` interface invokes this external program; the runner records its path, version, size and SHA-256.
- Python 3 and Rust/Cargo to build the existing source-bound rational point validator.

Do not provide a Magma passfile, license contents, tokens, or a process-environment dump as an artifact. The runner records executable identities and host/CPU fields without serializing the environment.

## Preflight and exact inventory

From a clean checkout containing merge `5dbe54dd468e914b5872d10d06e7c44d01b0a236`:

```bash
cargo build --locked --release --example koblitz_pdp_backend
python3 scripts/run_koblitz_external_magma.py --self-test
python3 scripts/run_koblitz_external_magma.py --dry-run-manifest > external-magma-plan.json
```

The dry-run manifest must report schema `koblitz_external_magma_dry_run_manifest.v1`, twenty tasks in seed-major order, Stage 13 source revision `67cf7f559cf6cf4c70fc16f22eb5dd68e4883d31`, and archive commit `5dbe54dd468e914b5872d10d06e7c44d01b0a236`. It also binds the Stage 13 verifier, matrix parser, meter, protocol, custody record and archived dependency tree to the pinned commit. Review the absolute paths and `relevant_checkout` state before execution.

## Full run

Use a new output directory and explicit executable paths:

```bash
python3 scripts/run_koblitz_external_magma.py \
  --magma /absolute/path/to/magma \
  --minisat /absolute/path/to/minisat \
  --backend /absolute/path/to/target/release/examples/koblitz_pdp_backend \
  --output /new/path/koblitz-external-magma-run
```

The frozen defaults are one worker, one requested Magma/F4 thread, GPU disabled by the input, 120-second watchdogs for primary F4, F4-plus-`SAT(G)` witness extraction and point validation, and at most 64 separately charged model attempts. A timeout, missing marker, proper basis without a validated rational point witness, source-valid nonlifting model, or model cap is inconclusive. A claimable full run also requires the pinned archive contracts and a clean relevant checkout; an uncommitted runner or altered dependency is recorded and cannot set the full-gate boolean.

For each task the runner:

1. Executes the byte-identical archived `instance.magma` through `process_meter.py`.
2. Parses the direct sparse F4 status, step degrees, basis size and internal timers.
3. For a proper basis, creates a deterministic second script from the exact Boolean source prefix. That script recomputes the same F4 basis and calls `SAT(G)`; its process receipt therefore charges both F4 recomputation and model extraction.
4. Rejects disagreement between the primary and repeated F4 terminal identities.
5. Sends every extracted assignment to a separately metered `koblitz_pdp_backend validate-model` process. That validator regenerates the algebraic source, checks the assignment, constructs rational factor-base lifts and requires a signed point sum equal to the public target.
6. Excludes source-valid nonlifting assignments and repeats, up to the recorded model cap. Duplicate or malformed assignments fail closed.
7. Rehashes every frozen input after execution and writes atomic task, progress and aggregate receipts.

Resume an interrupted run only with the identical tool identities and policy:

```bash
python3 scripts/run_koblitz_external_magma.py \
  --magma /absolute/path/to/magma \
  --minisat /absolute/path/to/minisat \
  --backend /absolute/path/to/target/release/examples/koblitz_pdp_backend \
  --output /new/path/koblitz-external-magma-run \
  --resume
```

Resume revalidates saved commands, raw output, process metrics, witness scripts, assignment files and validator reports. It refuses a changed leaf instead of silently skipping it.

An `--only-task seed-N/CELL` run is an operational smoke test. It cannot set the full-matrix gate. The aggregate summary reports separate and combined core-seconds, sequential wall time, and maximum process `ru_maxrss` for primary F4, F4-plus-`SAT(G)`, validators, and any recovered interrupted attempts, plus per-cell distributions. A synthetic fake-tool run requires explicit `--synthetic-test-mode` and is unconditionally gate-ineligible.

## External return packet

Return the complete output directory without editing it, together with:

- the saved dry-run manifest;
- the exact repository commit and clean/dirty status;
- reviewer identity or stable handle, affiliation, independence statement and conflicts;
- Magma license/access class without passfile or license contents;
- a source-pinned novelty review following `EXTERNAL_NOVELTY_REVIEW_TEMPLATE.json`;
- an explicit `CONCUR`, `QUALIFIED`, or `BREAKS` verdict and the narrowest supported claim.

Executable hashes and strict version probes bind the observed tools; they do not authenticate license entitlement or reviewer independence. A successful licensed-host matrix closes only the missing Magma execution/accounting portion. It does not by itself establish independent methodology, novelty, end-to-end index-calculus scaling, or a SOTA result.

Primary documentation: [Magma direct F4 and Boolean polynomial rings](https://magma.maths.usyd.edu.au/magma/handbook/text/1314) and the [Magma SAT/minisat interface](https://magma.maths.usyd.edu.au/magma/handbook/text/1318).

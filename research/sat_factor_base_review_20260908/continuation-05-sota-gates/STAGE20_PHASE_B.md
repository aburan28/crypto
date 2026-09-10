# Stage 20: blinded balanced PDP Phase B

This stage runs the 160 Phase-A explicit targets through native XOR, WDSat,
and CryptoMiniSat without exposing the Phase-A ground-truth ledger to the
solver driver.  It creates 480 solver outcomes.  Timeouts and `Unknown`
remain inconclusive.

The production run has not been started.  Build and review the committed
executables first:

```text
cargo build --release --bin koblitz_pdp_prepare --example koblitz_pdp_export --example koblitz_pdp_backend
```

That convenience command is currently **not** a charged production receipt.
The frozen protocol therefore keeps `full_cost_gate_passed` false.  The same
applies to the prebuilt CryptoMiniSat executable: its binary is hash-bound for
the panel, and its clean source is pinned at commit
`7ae1b4a74259cdce223a584281fb8f090bbd3eed`, but a fresh metered configure and
build receipt has not yet been implemented.  A complete solver panel remains
an internal engineering measurement until additive Rust and CryptoMiniSat
build receipts are frozen and bound; compilation cost must remain separate
from algorithmic solver time.

Build the frozen WDSat binary from the clean, pinned checkout.  This command
charges source copying, config installation, cleaning, compilation, and the
final binary copy separately.  The corrected archived n59 config is a frozen
overprovision candidate; the later export-only pass checks its exact capacity
before the first solver starts.

```text
python3 scripts/build_koblitz_phase_b_wdsat.py \
  --source /Volumes/SSD990/wdsat-binary-sat-review-20260909 \
  --output /Volumes/SSD990/koblitz-balanced-pdp-phase-b-wdsat-20260910
```

Create a sanitized solver root.  The staging command authenticates the raw
Phase-A seal and blind bundle, then copies only `blind/bundle.json` and a
path-free `binding.json`.  It does not copy the raw seal or ground-truth file.

```text
python3 scripts/run_koblitz_blind_pdp_phase_b.py stage \
  --phase-a-seal /Volumes/SSD990/koblitz-balanced-pdp-phase-a-20260910/prepared/seal.json \
  --blind-bundle /Volumes/SSD990/koblitz-balanced-pdp-phase-a-20260910/prepared/blind/bundle.json \
  --output /Volumes/SSD990/koblitz-balanced-pdp-phase-b-solver-input-20260910
```

Run the driver under one outer process meter.  The driver itself meters every
explicit export, isolated source reconstruction, solver, and point-witness
validation in a fresh process.  It exports and freezes all 160 sources, proves
the WDSat config covers their exact maxima, and only then begins the 480
solver attempts.  Use a new output root; interrupted runs are retained as
incomplete and are never resumed or overwritten.

```text
mkdir -p /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910

python3 scripts/process_meter.py \
  --cwd /Volumes/SSD990/crypto-koblitz-balanced-pdp-phase-b-20260910 \
  --timeout 86400 \
  --stdout /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910/driver.stdout \
  --stderr /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910/driver.stderr \
  --metrics /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910/driver.metrics.json \
  -- \
  /Library/Frameworks/Python.framework/Versions/3.13/bin/python3 \
  /Volumes/SSD990/crypto-koblitz-balanced-pdp-phase-b-20260910/scripts/run_koblitz_blind_pdp_phase_b.py run \
  --solver-root /Volumes/SSD990/koblitz-balanced-pdp-phase-b-solver-input-20260910 \
  --output /Volumes/SSD990/koblitz-balanced-pdp-phase-b-run-20260910 \
  --exporter /Volumes/SSD990/crypto-koblitz-balanced-pdp-phase-b-20260910/target/release/examples/koblitz_pdp_export \
  --backend /Volumes/SSD990/crypto-koblitz-balanced-pdp-phase-b-20260910/target/release/examples/koblitz_pdp_backend \
  --wdsat /Volumes/SSD990/koblitz-balanced-pdp-phase-b-wdsat-20260910/wdsat_solver \
  --wdsat-build-receipt /Volumes/SSD990/koblitz-balanced-pdp-phase-b-wdsat-20260910/receipt.json \
  --cryptominisat /Volumes/SSD990/cryptominisat-koblitz-20260909/build/cryptominisat5
```

Only after `run-seal.json` exists should the separate scorer receive the raw
seal and ground-truth ledger.  It independently replays the run inventory,
task count, 480 backend outcomes, resource totals, WDSat capacity, SAT witness
admission, and exact blind-id/target/cell join before emitting class labels.

```text
python3 scripts/score_koblitz_blind_pdp_phase_b.py \
  --solver-root /Volumes/SSD990/koblitz-balanced-pdp-phase-b-solver-input-20260910 \
  --run-root /Volumes/SSD990/koblitz-balanced-pdp-phase-b-run-20260910 \
  --outer-metrics /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910/driver.metrics.json \
  --phase-a-seal /Volumes/SSD990/koblitz-balanced-pdp-phase-a-20260910/prepared/seal.json \
  --oracle-ledger /Volumes/SSD990/koblitz-balanced-pdp-phase-a-20260910/prepared/sealed-oracle/oracle-ledger.json \
  --output /Volumes/SSD990/koblitz-balanced-pdp-phase-b-score-20260910
```

This is solver-view isolation: child arguments, environment, standard input,
working inputs, and task directories contain no ground-truth labels or paths.
It does not claim protection against a malicious executable that traverses
other host-readable paths.  That stronger claim requires a filesystem sandbox
or container exposing only the sanitized root, tools, and fresh output mount.

Even a complete result is a balanced public toy PDP measurement.  It does not
measure natural target prevalence, relation collection, linear algebra,
unknown-scalar discrete-log recovery, Pollard-rho crossover, external
reproduction, novelty, or Koblitz index-calculus SOTA.

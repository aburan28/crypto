# Stage 20: blinded balanced PDP Phase B

This stage runs the 160 Phase-A explicit targets through native XOR, WDSat,
and CryptoMiniSat without exposing the Phase-A ground-truth ledger to the
solver driver.  It creates 480 solver outcomes.  Timeouts and `Unknown`
remain inconclusive.

The production run has not been started. A normal development build is useful
for tests, but is not a charged production build receipt:

```text
cargo build --release --bin koblitz_pdp_prepare --example koblitz_pdp_export --example koblitz_pdp_backend
```

Use the measured builder for the Rust exporter/backend and CryptoMiniSat.
Each output directory must be new. The source checkout must be clean. The
builder archives the committed source, measures extraction, tool-version
probes, configuration, compilation and final executable copies, and retains
the raw process receipts. Cargo's ordinary hard-linked example outputs are
copied into dedicated single-link `bin/` artifacts before identity binding.

```text
python3 scripts/build_koblitz_phase_b_tools.py rust \
  --source /Volumes/SSD990/crypto-koblitz-balanced-pdp-phase-b-20260910 \
  --output /Volumes/SSD990/koblitz-balanced-pdp-phase-b-rust-build-20260910 \
  --jobs 2

python3 scripts/build_koblitz_phase_b_tools.py cryptominisat \
  --source /Volumes/SSD990/cryptominisat-koblitz-20260909 \
  --output /Volumes/SSD990/koblitz-balanced-pdp-phase-b-cms-build-20260910 \
  --jobs 2
```

The repository's root `Cargo.lock` is ignored. The build instead installs the
committed `stage-20-rust-build/Cargo.lock` snapshot in a metered step. Both
vendoring and compilation use `--locked --offline`, so an incompatible lock
stops the build. Vendoring reads the existing source cache; compilation uses
a new Cargo home, the recorded vendor configuration, and a new target
directory. Previously compiled objects are not reused. The receipt binds
the committed manifest, lock snapshot, library tree and example sources.

CryptoMiniSat is fixed at `7ae1b4a74259cdce223a584281fb8f090bbd3eed`.
The recipe separately archives the clean CaDiCaL and CaDiBaCk sources already
present in that checkout at `878a600ea92ce3a9173b069ab294429c851acf85` and
`acc8bb55b98e35afe0dc228662bf41d8341b8a9b`. CMake network fetching is disabled.
The eleven fixed test/utility submodules are recorded as unused exclusions;
the recipe disables testing and builds only `cryptominisat5`. Compiler paths,
tool versions, environment, CMake cache, source archives and binary hashes
remain in the artifact. Previously installed system libraries, including GMP,
are an explicit dependency of this build environment.

The original protocol stays byte-identical. The additive receipts are bound
by the execution plan, copied with their raw evidence into the terminal run
archive, and checked again by the separate scorer. Full production execution
and scoring require both receipts. A caller using `--allow-dirty` always gets
the `operational_smoke` evidence class, including on a clean checkout.

Build core-seconds include waited compiler/linker descendants and parallel
jobs. They are charged once per tool build, separately from instance export,
solving and validation; they are not repeated for every target. Summed stage
wall time is not whole-driver wall time, and core-seconds are not observed
single-core elapsed time. Largest-child high-water RSS is reported explicitly;
simultaneous aggregate build-tree peak RSS and measured single-core elapsed
time remain unavailable. The legacy meter field `single_core_seconds` is an
alias for CPU-seconds and is not promoted to a single-core benchmark.

For a whole-driver timing record, invoke each builder under
`scripts/process_meter.py` with separate stdout, stderr and metrics files.
That enclosing receipt includes the child stages, so its CPU time must not
be added to their sum. The builder receipts themselves exclude outer Python
overhead, prior source acquisition, and prior toolchain/system dependency
installation. These limitations retain `full_cost_gate_passed: false` even
after both build receipts are bound. Local builds also supply no remote Magma
process-memory, parallel-resource, licensing-cost, or external-review evidence.
Failed attempts retain their raw artifacts and cannot be resumed or overwritten.
The two enclosing build-driver receipts remain separately retained operator
evidence: the current Phase-B run and scorer do not authenticate or incorporate
them, so they cannot close the full-cost gate.

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
  --timeout 172800 \
  --stdout /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910/driver.stdout \
  --stderr /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910/driver.stderr \
  --metrics /Volumes/SSD990/koblitz-balanced-pdp-phase-b-outer-20260910/driver.metrics.json \
  -- \
  /Library/Frameworks/Python.framework/Versions/3.13/bin/python3 \
  /Volumes/SSD990/crypto-koblitz-balanced-pdp-phase-b-20260910/scripts/run_koblitz_blind_pdp_phase_b.py run \
  --solver-root /Volumes/SSD990/koblitz-balanced-pdp-phase-b-solver-input-20260910 \
  --output /Volumes/SSD990/koblitz-balanced-pdp-phase-b-run-20260910 \
  --exporter /Volumes/SSD990/koblitz-balanced-pdp-phase-b-rust-build-20260910/bin/koblitz_pdp_export \
  --backend /Volumes/SSD990/koblitz-balanced-pdp-phase-b-rust-build-20260910/bin/koblitz_pdp_backend \
  --rust-build-receipt /Volumes/SSD990/koblitz-balanced-pdp-phase-b-rust-build-20260910/receipt.json \
  --wdsat /Volumes/SSD990/koblitz-balanced-pdp-phase-b-wdsat-20260910/wdsat_solver \
  --wdsat-build-seal /Volumes/SSD990/koblitz-balanced-pdp-phase-b-wdsat-20260910/build-seal.json \
  --cryptominisat /Volumes/SSD990/koblitz-balanced-pdp-phase-b-cms-build-20260910/bin/cryptominisat5 \
  --cryptominisat-build-receipt /Volumes/SSD990/koblitz-balanced-pdp-phase-b-cms-build-20260910/receipt.json
```

The 48-hour outer watchdog covers the frozen sequential worst-case envelope:
160 exports, 160 isolated-source checks, 480 solver processes, and up to 320
conditional external witness validators, each with a 120-second child
watchdog, plus driver overhead. The outer receipt encloses those child costs
and therefore is retained separately rather than added to their totals.

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

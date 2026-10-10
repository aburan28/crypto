# Chained-S3 Linux worker startup gate

The offline `x86_64-unknown-linux-musl` build at source commit
`7eba91dd46661a2e035c6fe1f5eb15a2e09ba7fa` completed and produced a
statically linked ELF. `build-receipt.json` records its exact command,
compiler versions, executable SHA-256 and full `build.log` hash. The binary is
an external build artifact, not a committed repository file.

The no-network Docker smoke used a 256 MiB hard memory cap, zero swap, one
CPU, a 60-second wall cap, `--max-trials 0`, and deliberately invalid JSON for
the panel and public fixtures. The worker reached panel validation after its
cgroup and compiled-source checks. It rejected the invalid panel with exit 1;
the outer status is the expected `PRODUCER_FAILURE_worker_exit`. No events,
relation row, SAT search, rank stage, target scalar or total-runtime claim was
produced. This verifies launch and guard wiring only; it does not validate
larger-base import capacity or a natural N83 relation.

The archived `config.json`, `worker-config.json`, `outer.json`, source
attestation, input files and process logs bind the outcome. Run
`python3 research/koblitz_n83_factor_base_sweep_20261008/verification/chain-search-worker-startup/verify.py`
from the repository to check their hashes against the committed source. When
the local ELF still exists, that verifier also checks its bytes against the
build receipt. The original executable can be rebuilt from the recorded
offline command. A retained-base construction or search still requires a new
compute allocation and its own fresh output directory.

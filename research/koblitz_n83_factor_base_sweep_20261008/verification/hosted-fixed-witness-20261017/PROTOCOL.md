# Hosted N83 fixed-witness SAT-access gate

This single, bounded control tests whether the current chained-S3 solver can
recover the independently group- and S3-verified, ordered five-point witness
from the retained K1182 small domain and K2048 residual suffix. It fixes all
five summand x-coordinates and all three affine intermediate x-coordinates.
It is a solver and encoding access check, not a natural relation-yield trial.

The committed `domains.json` from the earlier control has SHA-256
`63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb`.
The workflow checks that digest before execution, checks all source and input
hashes again before solving, and checks out its dispatch commit by exact SHA.
The pinned Rust image is used on an ARM64 GitHub runner. The solve runs with
one CPU, a 4 GiB memory limit, no swap allocation, a 5,000-conflict solver cap,
and a 1,200-second external wall cap. Fetch and build have separate bounded
stages. The three phase markers in `solve.stderr` identify circuit completion,
domain completion, and fixed-unit installation even if the wall cap fires.

The artifact retains source and input digests, runner/image identity, build
logs, worker binary/digest, raw solve streams, exact exit code, wall bounds,
container inspection and cleanup. A successful `SAT_verified_witness_x`
requires `model_verified=true`, a checked x-coordinate match and original
S3-system replay by the worker. An `UNKNOWN_conflict_cap` or external exit 124
is censored. A reported `UNSAT_fixed_witness_encoding` contradicts the prior
algebra replay and requires a source/encoding investigation. An invalid SAT
model is a producer failure. None of these outcomes alone ranks factor bases.

The dispatch may consume at most 40 hosted runner minutes; the search cap is
20 minutes, and the build watchdog is 10 minutes. The later retained receipt
must separately report measured setup and solve walls, conservative charged
capacity, cleanup and audit. Source, build and search times cannot be silently
excluded from a complete IC comparison.

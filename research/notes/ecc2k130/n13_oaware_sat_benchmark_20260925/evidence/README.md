# Frozen source and raw evidence

The single base CNF is the merged #781 file at
`../rotated_s3_o_branch_20260925/evidence/producer/n13-m5/base.cnf` from
this study directory, 17,630,779 bytes, SHA-256
`ac6fb610d38023c392308826c88dd5b0c2bbe1347057bd3f29c6a3b9e05d468b`.
`FROZEN.json` pins it and the schema, independent point corpus, source,
interpreter/tool and three solver binary hashes. `INPUT.json` pins all
32 full-point Q+T assumptions and independent oracle labels. The source
commit used for outcomes is `cea74c8e0133b972d5cb0e2eff7de2c73feeb308`.

`smoke/` contains six tiny SAT/UNSAT child receipts and raw output.
`panel/` contains 96 JSON child receipts and stdout/stderr bytes; each
receipt names the one temporary derived `query.cnf` file, its full-file
SHA-256 and byte count, positive unit assumption and exact cold command.
The temporary 17.6 MB query file was deleted after each target's three
children, so no 96-copy base is committed. `ci_replay.py` reconstructs
all 32 query bytes from the one merged base and rechecks each solver
status and signed exact-point witness. `analysis.json` provides the
preregistered fixed-order Q and full-portfolio costs.

`preoutcome_setup_failure.txt` records the Python 3.9 import failure;
`preoutcome_smoke_invocation_failure.*` records the wrong 3.12 shim;
`failed_smoke_sandbox/` and `failed_smoke_sandbox.stderr.txt` retain the
first exact-version tiny smoke aborted by sandbox-denied psutil sysctl.
These are setup failures, not corpus results. `smoke-run.*` and
`panel-run.*` retain top-level stdout/stderr. `MANIFEST.json` lists
SHA-256 and byte count for every raw receipt/output and analysis artifact.

Replay from repository root, without a solver process:

```sh
python3.12 research/notes/ecc2k130/n13_oaware_sat_benchmark_20260925/ci_replay.py \
  --smoke research/notes/ecc2k130/n13_oaware_sat_benchmark_20260925/evidence/smoke \
  --panel research/notes/ecc2k130/n13_oaware_sat_benchmark_20260925/evidence/panel
```

The replay validates SAT models by independent point-fibre/group law and
compares UNSAT statuses with the exhaustive toy point oracle. It does not
check a separate UNSAT proof certificate. All timing is a toy solver-stage
measurement; no matched rho or full ECDLP cost is claimed.

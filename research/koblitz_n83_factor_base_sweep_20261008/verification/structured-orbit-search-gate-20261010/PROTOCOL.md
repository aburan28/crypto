# N83 structured-orbit import and first search gate

**Status: executed once; see [the result](RESULT.md).** The prior local pilot ledger
charged all 7,200 of the previously authorized seconds. The user extended the
limit for this gate by **1,800 additional local seconds**:
600 seconds for the capped offline Linux build and 1,200 seconds for the two
experiment stages and cleanup. The authorization receipt must name this gate,
the exact fresh build and experiment directories, the earlier 7,200-second
charge, the new 1,800-second allowance and a reference to that approval. The
build and search supervisors reject a missing or insufficient receipt and
charge their full 600 and 1,200 seconds after a launch, even if a stage ends early.

## Frozen input and workload

Use the downloaded, byte-verified object from
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v3-structured-orbit/a0/objects/aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389.jsonl.gz`.
Its compressed SHA256 is
`c5de0b26e84ad80bfb7be6a214b89c9682e356effb6079c7abf28cd1151ffb03`.
The manifest is the committed `structured-d12-orbit-20261010/manifest.json`;
the public target is `a=0, fixture=0` in `pilot-01/probe-corpus.json`.
`probe-validation.json` is mounted separately and opened by the worker only
after a candidate scalar exists. The worker is
`examples/koblitz_n83_structured_orbit_worker.rs` compiled as a release Linux
ELF from a clean source commit. The adjacent `build-linux.sh` freezes the
offline arm64 build command, caps compilation at 450 seconds with the rest of
its 600-second charge reserved for cleanup, and writes a binary/source/image receipt outside
the checkout. Record its SHA256, source commit, source-file hashes, local
Docker image ID, input hashes and exact command before launch. The build
script was executed once for the worker; the frozen build receipt is retained
under `evidence-01/build/`.

`structured_orbit_search_supervisor.py` launches two fresh no-network Linux
containers with one CPU, 4 GiB memory, zero swap, 64 PIDs, dropped capabilities
and a read-only filesystem. The first has a 240-second outer wall cap and runs
only `verify`: it checks the archived source projection, all 332,166 labeled
point rows and the reconstructed 2,001-column generic base. A failed or capped
import prevents the next stage. The second has a 720-second outer wall cap and
runs one five-summand `ChainedS3FrobeniusOrbit` natural trial. The worker's own
deadline is 710 seconds; its SAT call has at most 10,000 conflicts and 32
models, 150,000 variables and 12,000,000 domain clauses. Container cleanup is
checked after every launch; at least 120 seconds of the 1,200-second charge
are reserved for cleanup and receipts. A timeout, memory cap, malformed receipt, source
change or cleanup failure cannot be promoted to relation evidence.

The output is a fresh directory outside the checkout. It keeps `config.json`,
stage stdout/stderr, `outer.json`, `budget.json`, and the worker's progress,
attempt, relation, rank and cgroup peak-memory records if the search stage
reached them. The search supervisor defaults to the same immutable Linux image
ID pinned by the build script. A prior local `python:3.11-slim` tag resolved to
a different image and would fail the frozen build/search identity check before
archive import.

The outer receipt validates the worker's counts and retains the status
`UNKNOWN_solver_cap`, `UNKNOWN_trial_cap`, `UNKNOWN_wall_cap`,
`UNKNOWN_resource_cap` or `PRODUCER_FAILURE_*`. One trial cannot determine the
2,001 relation columns and the driver disables direct-relation bypass; a
purported verified target in this gate is treated as a receipt anomaly for
independent investigation. A natural relation, if found, still requires
independent group replay. This gate does not establish full precomputation,
an individual log stage, matched cold cost or the best factor base.

## Authorized launch procedure used

The approved launch used a local `authorization.json` outside the checkout
with schema `n83.additional-local-budget/v1`, study
`koblitz_n83_factor_base_sweep_20261008`, purpose
`structured_orbit_import_and_one_m5_trial`,
`prior_pilot_charged_seconds: "7200"`,
`additional_authorized_seconds: "1800"`, `build_dir` and `output_dir` set to
the exact fresh absolute build and experiment paths, and
`user_approval_reference` identifying that approval.
Do not commit that local authorization reference. From the clean committed
checkout, build the Linux worker in a fresh directory, then use that binary
and the downloaded object:

```sh
sh research/koblitz_n83_factor_base_sweep_20261008/verification/structured-orbit-search-gate-20261010/build-linux.sh \
  /absolute/path/to/fresh-build /absolute/path/to/authorization.json
python3 research/koblitz_n83_factor_base_sweep_20261008/structured_orbit_search_supervisor.py \
  --object /absolute/path/to/aea278aff2d79a592d5164d7437ddde8f8e9377fa481058dfc16afa176a31389.jsonl.gz \
  --binary /absolute/path/to/fresh-build/worker \
  --authorization /absolute/path/to/authorization.json \
  --output /absolute/path/to/fresh-output
```

The commands ran once against the frozen source and object. The resulting
build, import, search, cleanup and budget receipts are under `evidence-01/`.
The supervisor's synthetic controls also validate the missing-budget guard,
archive receipt rejection, resource command and prevention of search after a
failed import; the real bounded result is reported separately in `RESULT.md`.

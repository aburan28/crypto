# Frozen-admission replay across later source changes

Status: preregistered repair, round 2. This protocol and its checker are
frozen before the new current-tree contract result is read. It does not
rerun a PDP, admit a solver, or change the 2026-09-25 admission verdict.
Round 1 (below, "Trigger and historical claim" through "Bounded run and
decision") covers the builder-only drift and is unchanged; round 2 ("A
second reference file drifts: from one snapshot to a snapshotted set")
covers a second file drifting two days later and generalizes the
mechanism from a single snapshot to a snapshot per drift-exposed file.

## Trigger and historical claim

The original `rotated_solver_admission_20260925/INPUT.json` records commit
`8b640f31195238242be8940c928454ad2c9f2ce6` and ten exact reference
SHA-256s. On current main at `8b1e2d3734c1f30c14ce253b33094f894ab5005d`,
the old `ci_replay.py` fails its reference-file hash assertion: nine reference
files still match and `src/cryptanalysis/koblitz_groebner.rs` changed from
`bd95ee98ad095f93aec5a2b79e7f4840dcf4473db174418aef9b77ebae79c015`
to `429be4bcfb520e12481d1f95a7c54862cab3f691ab585adbb60be1ba7d033518`.
The later changes add checked layout arithmetic and explicit unsupported-input
reporting. This is a replay-environment failure, not a changed historical
experiment result.

The original `INPUT.json`, `PROTOCOL.md`, `FROZEN.json`, `audit.py`, `run.py`,
`ci_replay.py`, and both raw evidence directories are immutable. A byte-for-byte
copy of the historical builder from that recorded commit was stored, in this
round, as `historical_koblitz_groebner.rs.gz` (superseded in round 2 below by
`historical_snapshots/koblitz_builder.rs.gz`, identical bytes, moved to make
room for the other four). The checker verifies its *decompressed*
SHA-256 against the original frozen `INPUT.json`; it also verifies all nine
other reference files against their original frozen hashes. It constructs an
isolated temporary tree with those ten pinned reference bytes and the five
frozen runner files, then invokes the **original** `ci_replay.py` there both
without and with the original final receipt. The original replay must verify
the existing blocker result unchanged. GitHub CI needs no checkout of an
unadvertised old commit; source provenance is the recorded commit plus the
verified byte-for-byte snapshot. The fixture is compressed only to keep this
repair small, not to change the bytes the old audit sees.

Failure of any historical hash, missing file, subprocess, receipt, or old
assertion is a failure. No new result can substitute for the old receipt.

## Separate current-tree interface contract

The current source is checked *as current source*, with its digest reported,
not forced to equal the old digest. This bounded contract concerns the
generic `build_decomposition_system`, the actual
`polynomial_reuse::build_decomposition_system_reusing`/`DecompositionTemplate::build`
path used by `groebner_decompose`, and that frontend's unsupported status:

1. `MAX_VARS` is 64; both the direct builder and cached template still use
   one `basis` for every summand. The reuse wrapper must instantiate that
   template or fall back to the checked direct builder.
2. Both builders reject `m < 2`, use checked multiplication and addition for
   the layout, and return `None` when `n_vars > MAX_VARS`. The template also
   rejects unsupported field widths. For the frozen toy arities, `n13,m5,d2`
   has 49 raw bits and fits; `n19,m6,d2` has 88 and does not.
3. The Gröbner frontend maps an unsupported layout or infinity target to
   `exhausted=true, unsupported=true`, rather than a refutation. The existing
   Rust test `pdp_admission_unsupported_is_not_a_refutation` and a new
   n19/m6/d2 width-only regression must pass. That regression uses a
   deliberately inert factor-base sentinel because rejection precedes every
   point-domain or witness operation; it does not test a rotated PDP.

Passing this narrow contract does **not** conclude that no new complete
rotated exporter exists elsewhere, nor does it certify O-aware PDP semantics
for this generic builder. A future code change that alters one of these
properties requires an explicit new contract assessment; the historical
replay still stands. The current source digests are diagnostics, not new
frozen input hashes.

## Bounded run and decision

After preregistration passes, run the historical replay and current contract
once. Both are deterministic read-only checks. No solver process, benchmark,
or binary inventory runs. Archive stdout/exit status, current source digests,
and the final decision in this PR. The path-filtered existing admission
workflow now runs the new freeze, historical and current checks, fail-closed
Python controls, and both targeted Rust regressions. Its bytes are part of
this new freeze. A shallow checkout is sufficient because the 52-KiB
historical source fixture is committed in this PR, decompressed under a size
cap, and checked against the old SHA-256. CI repeats the two checks and the
targeted Rust tests. Success means only `HISTORICAL_REPLAY_PASS` and
`CURRENT_GENERIC_CONTRACT_PASS`; a failure remains visible and is not relabeled
as a PDP result. There is no performance ratio or ECC2K-130 transfer claim.

## A second reference file drifts: from one snapshot to a snapshotted set

Two days after round 1 landed, `frozen-source-and-archive` (the original,
untouched `ci_replay.py`, run directly by CI on `main`'s tip against the
path-filtered workflow) failed again on `main`. Reproducing it directly
(`python3 research/notes/ecc2k130/rotated_solver_admission_20260925/ci_replay.py`
on `main` at `6c0324c9`) shows **two** of the ten reference files no longer
match the original `INPUT.json`, not one:

| reference key | path | frozen sha256 | current sha256 on `6c0324c9` |
|---|---|---|---|
| `koblitz_builder` | `src/cryptanalysis/koblitz_groebner.rs` | `bd95ee98…79c015` | `1590fd29…881a63c0` (drifted again, past round 1's `429be4bc…` snapshot) |
| `sat_example` | `examples/koblitz_s5_sat_instance.rs` | `c2bc8b05…2115f09` | `4ce1a95c…d550c8eaf` |

Both are ordinary files under `src/cryptanalysis/` and `examples/` that
unrelated cryptanalysis work keeps editing; nothing about either file is
specific to this admission thread. Waiting for each such file to break the
replay in turn, and hand-repairing it, is not a plan — it is round 1 repeated
indefinitely. Round 2 generalizes the round-1 mechanism instead of repeating
it: **every** reference-file key whose path lives under `src/` or `examples/`
(`binary_semaev`, `koblitz_builder`, `sat_example`, `solver_adapter`,
`wdsat_adapter` — five keys, confirmed by inspecting all ten paths in
`INPUT.json`) gets its own committed gzip snapshot at `BASE_COMMIT`, verified
byte-for-byte against `INPUT.json`'s original hash exactly as round 1 verified
the builder's. The remaining five reference files (`corpus_archive`,
`export_gate_result`, `fibre_theorem_result`, `n13_s3_result`,
`n19_s3_result`) are evidence artifacts inside this thread's own
`research/notes/ecc2k130/rotated_*_20260925/` directories, which nothing
outside that thread edits; they stay verified against the live tree, both
because that is cheaper (the archive alone is 2 MB) and because a live check
that keeps passing is itself evidence that no unrelated drift has reached
them.

Every one of the five historical snapshot bytes was checked against
`INPUT.json`'s recorded hash for that key before being committed — the base
commit's tree is the sole source, never `main`'s current tree — and the
gzip round-trips back to those exact bytes. `historical_replay` now reads
each `SNAPSHOT_KEYS` member from `historical_snapshots/<key>.rs.gz` and every
other reference key from the live tree, unchanged from round 1's per-file
branch. The size cap loosens from round 1's 1 MiB (sized for one file) to
4 MiB total headroom across the largest of the five (285,865 bytes,
`sat_example`), still far above any of them and still there only to bound a
truncated or corrupted stream, not to license a larger fixture.

This round changes no verdict, no protocol claim, and no original artifact;
it only widens which reference bytes are pinned by a snapshot instead of a
live check, because the round-1 assumption — that only the builder would
drift — held for two days and no longer does.

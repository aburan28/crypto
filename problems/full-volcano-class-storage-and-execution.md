# Full large-volcano class storage and execution

- **Status:** open
- **Kind:** resource limit; missing qualified graph-wide cryptanalysis execution;
  Conductor access/tool limit
- **Affected work:** [fixed five-class protocol](../research/full_volcano_classes_20261016/PROTOCOL.md)
- **Next action owner:** current implementation session for code and controls;
  full-run infrastructure owner unassigned
- **Last updated:** 2026-10-16

## Problem

The complete five-class request is not established because the remaining full
graphs need more durable random-access storage than this workspace provides,
and the graph-wide multi-method cryptanalysis workload has not been executed.

The exact five public input models are retained in the study's `inputs/`;
their original historical case labels are filenames, not curve identities.
Together their expected classes contain 7,833,658,967 vertices and 515 conductor
strata. The requested output includes every vertex, all rational kernels for
each selected finite prime-degree set, all strata and explicit maps, plus
verified graph-wide DLP results, operation distributions, timing and IC yields.
Those are separate acceptance conditions; merging tooling cannot discharge them.

## Evidence and uncertainty

- **Observed:** the workspace filesystem is 32 GiB; approximately 10 GiB was
  free after preserving the complete historical B graph and its compressed
  copy, checking out the repository and installing the native toolchain.
  The separate temporary mount is 8.8 GiB. The session has a 16 GiB memory
  limit and four-CPU quota. Dependency installation reached ENOSPC; removing
  only the newly downloaded unused Sage image recovered space. Dataset files
  and remote Drive copies were preserved.
- **Derived:** only three u64 vertex fields (a,b,j), with no indexes or edges,
  require 188,007,815,208 bytes across all five classes. The new SQLite format
  additionally stores kernels, scalings, indexes and journals; the 24-byte
  figure is a fixed-layout baseline, not a measured SQLite storage estimate or
  an information-theoretic bound. Compressed/generative encodings may be smaller.
- **Observed historical evidence:** B's native manifest records all 129,140,163
  vertices, 18 strata, 258,280,326 directed degree-3 kernel maps, 9,298,091,736
  raw bytes and 7,904,902,167 compressed bytes. Its frozen audit and Drive
  reassembly manifests are copied into the study's `evidence/historical-B/`.
  Historical runtime measurements were made outside this repository's
  ecbench/isolation protocol and are not newly qualified performance results.
- **Tool gap:** `conductor check --summary ... --scope path:...` returned
  exit 127 (`conductor: command not found`); the base tree has no
  `.conductor/WORKFLOW.md`. No reservation or task receipt exists.
- **Missing capability:** a provisioned durable volume/object-store-backed
  execution system large enough for full enumeration, replay and benchmark
  records; completed matched ecbench implementations/workloads for the full
  requested comparison; an independent implementation replay of the new
  generic database backend.
- **What this does not establish:** mathematical impossibility, absence of
  other curves or isogenies, or a bound on achievable streaming/compressed
  representations. A different external-memory backend could change storage
  needs. Current resource limits do not cancel the user's requested workload.

## What would unblock it

- [ ] Retain complete vertex/map outputs for each fixed case on durable
  storage, with checksummed manifests and accessible download locations.
- [ ] Exhaust the finite degree set at every vertex, replay every map and
  verify exact stratum populations; preserve failures and interruptions.
- [ ] Run every requested DLP arm on the specified vertex/target population,
  with matched inputs, charged setup/verification, complete operation counts,
  end-to-end timing and IC relation-yield records under ecbench admission.
- [ ] Report min/max/mean with the actual denominators and missing outcomes,
  plus independent replay; no sampled population is labeled graph-wide.

## Next discriminating step

Use the native per-case `plan` receipts to provision a durable random-access
volume. Generate and replay one full remaining class first, retaining its
exact degree set and stratum counts, then continue the other fixed cases.
The existing small exhaustive controls verify correctness but do not complete
that full-run step. A prefix is retained as partial. Recheck on receipt of an
existing volume/worker or a validated streaming backend.

## Resolution and history

2026-10-16: protocol and native CLI implementation started. Full A/C/D/E
generation and the requested graph-wide benchmark campaign remain open.

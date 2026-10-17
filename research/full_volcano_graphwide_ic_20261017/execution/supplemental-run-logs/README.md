# Supplemental logs, 17 October 2026

The new auxiliary run logs, audit stdout and duplicate per-job specs are in
`runs.tar.xz`. Extract with `tar -xJf runs.tar.xz`; check the archive with
`sha256sum -c runs.tar.xz.sha256`. Repacking the extracted tree with the
AGENTS.md §13a recipe and `--mtime='2026-10-17 00:00:00Z'` produces identical
bytes (verified locally).

The sealed five-file session capsules and hashed child stderr remain loose in
`../campaigns/`, along with standard/exhaustion replay receipts. The lab browser,
native verifier and checkpoint manifest name those individual paths; this is
the explicit §13a path-resolution exception. Existing committed control, pilot
and first-12 capsules are unchanged. Only the subsequent 87 survey capsules are
added under `survey-vertices-13-99`, so the browser does not index them twice.

The complete stopped 99-vertex campaign, including ledger, full queue, exact
binaries and source, is separately delivered on Drive; see
[`../delivery-99.json`](../delivery-99.json). The optimized control/pilot have
separate capsules and a separate native summary. No debug/release pooling is
performed.

The repaired control and pilot have four new byte-preserved loose capsules
under `../campaigns/repaired-control/` and `../campaigns/repaired-pilot/`.
The native verifier and browser resolve those paths, so the same exception
applies. Their actual native counts are in `../repaired-qa-summary-v2/`; they
are not pooled with earlier repeated workloads. C's standalone diagnostic
raw tree is separately archived at `../finite-field-k1-evidence/runs.tar.xz`.

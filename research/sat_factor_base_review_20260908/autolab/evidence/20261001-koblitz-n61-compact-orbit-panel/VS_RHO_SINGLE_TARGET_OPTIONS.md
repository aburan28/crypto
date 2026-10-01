# n=61 compact orbit: can `vs_rho` be filled honestly? (2026-10-01)

The `vs_rho` stage of `docs/ic/boundary_targets.json` requires one target per
workload, a target-online wall interval for each arm, an exclusive five-phase IC
split that sums to the IC online total, a rho policy record, matched resource
envelopes, replay certificates and IC1 candidate/workload/run identities. The
n=61 batch panels (PR #830, PR #1134 and the autolab run in this directory)
fill none of the single-target fields. The claim drafts leave them absent, so
claim-check fails closed. This note lists what the current producers can supply
and what would need a decision. **No option is chosen here.**

## What the producers already emit

| Requirement | Compact orbit (`koblitz_orbit_dlp_fast`) | Rho |
|---|---|---|
| One target | Yes. A one-line scalars or `[x,y]` file works; setup (base, index, rank, LA) is still paid in full | KS v2 (`koblitz_rho_batch_ks_v2_n61`) and strong batch (`koblitz_rho_batch_ks_strong`) both accept `<fixtures>=1` at n=61 |
| Online start/stop | Per target: `target_ms` runs from the query hash to the end of the `[d]G = Q` check. Fixture construction is excluded (`target_generation_ms_excluded`) | Strong batch: per-fixture `total_ms` starts **before** Q is built from the known scalar, so it charges fixture construction to rho. KS v2: only process-level `setup_ms` / `in_process_ms` |
| Five exclusive IC phases | **Four** timers: `target_query_ms`, `target_pdp_and_relation_check_ms` (PDP and relation check fused), `target_descent_ms`, `target_recovery_check_ms`. Their sum is `target_phase_sum_ms`, a few µs short of `target_ms` | n/a |
| Rho policy | n/a | Strong batch: `lanes`, `dp_bits`, `table_entries`, `table_payload_lower_bound_bytes`, rung; single worker. KS v2: `dp_bits`, `table_payload_lower_bound_bytes` |
| Admissible reference | n/a | `rho_reference_minimum` (since 2026-10-01) admits `koblitz_rho_fixture … strong` (rejects n=61: `matches!(n, …, 53 \| 59)`), `koblitz_rho_batch_ks_strong` rung 3, and `koblitz_rho_batch_ks_v3` normal-basis. **KS v2 is not listed**, so it would need a same-target per-step cost no worse than strong |
| Scalar verified | `group_verified`, `recovered_matches_published` | `verified`, recovered == planted |
| Replay certificate | `independent_replay.py` receipts (pure Python), hashable | `growing_n_n61_L65536_20260930_rho_replay.py` receipts, hashable |
| IC1 identities | None. The harness would have to derive candidate/workload manifests | None |
| Resource envelope | Measured peak RSS ~2.5 GiB at K=600 | Measured peak RSS in MiB. Equal *caps* (threads=1, a common memory cap) would match; equal *usage* never will |

## The gaps that are not just plumbing

1. **The fused PDP + relation-check timer.** The schema wants
   `T_target_PDP_ms` and `T_target_relation_check_ms` separately, and the IC
   phase values must sum to `ic_online_wall_ms`. `extract()` interleaves
   candidate generation and the relation check inside one probe loop, so no
   honest split exists in the current output. Writing 0 for one phase, or
   splitting the fused value by a ratio, would be fabrication.
2. **The rho comparator.** The campaign's frozen comparator (KS v2) is not on
   the admissible list. The only listed single-target reference that runs at
   n=61 is `koblitz_rho_batch_ks_strong` rung 3 at L=1, and its per-fixture
   interval includes building Q from the known scalar.
3. **Identities and envelope.** The IC1 candidate/workload/run identity scheme
   and the definition of "same resource envelope" (caps or measured usage) are
   policy choices for the ledger owners, not for this harness.

## Options (for a producer/schema decision)

**A. Instrument the compact producer.** Add a split timer inside `extract()`:
accumulate PDP-candidate time and relation-check time per probe and emit
both, plus an explicit online start/stop event pair. Pair it with
`koblitz_rho_batch_ks_strong` rung 3 at L=1, changed to start `total_ms` after Q
is built, or to emit a separate `walk_ms`. This is a producer change to two
frozen examples (new source hashes) and needs a fresh one-target panel. It is
the only route to a schema-complete PASS without a schema change.

**B. Amend the schema for fused phases.** Allow an IC phase list that declares
which schema phases a timer covers (for example
`T_target_PDP_and_relation_check_ms` covering both), as long as the declared
phases still sum exactly to the online total. This needs no IC producer change,
but it changes the fail-closed contract. Ledger owners decide.

**C. Add a single-target mode to the batch harness.** Add a `--targets 1`
path that runs W one-target workloads, each with fresh rho state, and re-pays IC
setup outside the online interval. This is cheap to build on top of
`launch-panel`, but it still needs A (or B) for the phase split and a rho with
a clean online start. On its own it cannot pass.

**D. Keep `vs_rho` fail-closed for this campaign.** Report the batch panels as
`multi_target_batch_diagnostic` under `end_to_end_dlp`, as now, and record that
the operative reference (strong rho R3) gives IC/rho = 2.297 in retired
instructions at n=61 L=1,024
(`research/notes/index-calculus/RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md`).

Under A or C, any one-target n=61 comparison should be read against strong
rho, not KS v2: at L=1 the shared-setup amortisation that drives the batch
lead disappears, and IC pays its full K²·n index build for one target.

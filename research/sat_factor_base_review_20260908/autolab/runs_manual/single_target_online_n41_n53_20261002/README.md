# Single-target online `vs_rho` records — n=41 and n=53 (2026-10-02)

Committed evidence for the first two rungs of the primary single-target online
ladder on the Koblitz regime (`koblitz.vs_rho.n41_charged`, `koblitz.vs_rho.n53_charged`).

## What is claimed

- Exactly one previously unseen target per run; the IC arm (`koblitz_rank_fixture`)
  and the rho arm (`koblitz_rho_fixture`) solve the **identical frozen public point**
  under the same 16 GiB resource cap.
- Online clocks start **after** reusable IC preparation (factor base, support
  index, shared factor logs) and at rho's first target-dependent walk step;
  process launch, input loading, fixture generation, and reusable IC setup are
  excluded from both online intervals.
- n=41: three paired runs, online speedup **17.7x–18.8x**; both arms verified on
  every run; standalone Python GF(2^41) scalar replay PASS on two Python builds.
- n=53: two paired runs, online speedup **51.0x–55.7x**; both arms verified on
  every run; standalone Python GF(2^53) scalar replay PASS on both runs
  (`n53/replay_run_gf2n.py` derives curve, modulus, and scalars from the run logs).

## Explicit non-claims

Public synthetic fixtures only. Constant-factor wins only — no asymptotic
sub-rho claim, no key recovery, no external/private targets, no deployed-curve
security impact. The n=53 rung precompute materializes the explicit pair table
(peak child RSS 9.34 GB, under the 16 GiB cap, excluded setup, logged); the next
rung (n=61) must drop it for the compact-orbit domain. Multi-target amortized
results remain secondary and are not promoted by these records.

## Files

- `results.json` — summary record for both rungs.
- `n41/<run-id>_claim_draft.json` — per-run paired claim drafts (online phase
  costs, paired target, verified scalars, whole-process wall for reference).
- `n41/<run-id>_validation.json` — standalone Python GF(2^41) scalar replays.
- `n41/rejected_pairing_audit_20261001T020900Z-d7138bdf44.json` — the earlier
  schema-only PASS whose arms used different public points (`PAIRING_REJECTED`);
  retained as the audit that motivated the paired contract.
- `n53/<run-id>_claim_draft.json`, `n53/<run-id>_validation.json` — same contract
  at n=53 on the retained base `d859319…`.
- `<run-id>_state.json`, `*_resource.json` — autolab run states and resource
  receipts.

Full producer logs and receipts remain under the (gitignored)
`research/sat_factor_base_review_20260908/autolab/runs/<run-id>/` on the
measurement host; the censored failed attempt
`20261002T213515Z-14ab45f3ac` (precompute rank starved at 188/221 with the
eta 1/16 pointwise knobs) is retained there as well.

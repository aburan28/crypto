# Round 0004 pre-registration

Parent: round-0003b selection combined_batch8; source `/home/ubuntu/crypto/research/ic_candidate_tournament_20260915/runs/round-0003b/source_candidates/combined_batch8/source`. No parent holdout
data was used to select these mechanisms. Execution requires that same candidate
to pass the parent's full confirmation and replay.

Fresh seed 2026091604. Five challengers: fast factor-base lifting, folded pair
table, both, both with dense scalar algebra, both with batch4. Base support,
decomposition size, target count (one), all verification and full costs remain
fixed. Mathematical coverage and weak K-instruction floor stay unchanged.
Class: engineering. Same >=20% instruction and native-time promotion gates,
paired 95% limits, per-cell limits, confirmation and replay; same explicit
candidate/rho parity criterion. Budget 1800 paired jobs.

Builds are restricted to CPUs 0–1 while the parent measures on CPU 7. Round4
measurements start after the parent finishes, with every arm pinned to CPU 1.
The parent's full replay supplies a timing check after these builds complete.
Each round's runtime claims use its own fresh matched reference/candidate runs.

Preflight evidence: `preflight-next-tests.log` compares fast lifting against
general arithmetic (including point order) and folded decomposition against an
explicit full table on all five curve cells. These are correctness tests, not
performance claims. Frozen E2E oracle checks still cover every measured input.

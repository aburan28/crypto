#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.3pt)

= Binary regime A: the first sub-$2^(2ℓ)$ decomposition oracle at $ℓ = 8$

2026-10-15. `research/binary_mitm_oracle_l8_20261015/`; code `examples/binary_mitm_oracle_l8.rs`.
*Class: engineering* (per-target oracle work falls by $2^ℓ$; the $Θ(2^n)$ relation-collection total and the $ρ$ frontier are unchanged).

== The open row this closes

The boundary ledger's binary (non-Koblitz) `decomposition` row and agent priority #4 asked for the
first sub-$2^(2ℓ)$ oracle at fixed $ℓ = 8$ on $≤ 64$ targets, with median cost at most half of
pairs-and-solve, zero disagreements against the group check, and the measurement-schema fields
(including an explicit FFD status). The SAT/SMT survey's own falsification note
(`RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md`, §6) recorded the sub-$2^(2ℓ)$ question as "open and
un-attempted by anything in §6". This run is that oracle's first landing.

== Instance

- Field $F_(2^25)$, sparse Rabin-checked modulus (`find_irreducible_sparse(25)`). $n$ odd so the
  half-trace lifts $y$; $3ℓ = 24 ≤ 25$ keeps a genuine found/refuted mix.
- Curve $y^2 + x y = x^3 + a x^2 + b$ with $a = z ∉ F_2$ (genuinely non-Koblitz) and $b = 1$.
- Factor base: all curve points over the low-order subspace $⟨1, z, …, z^7⟩$ — 256 points
  (128 abscissae, both signs).
- Reference arm: `SubspaceOracle::low_order(8, 1)` — Gaudry pairs-and-solve over the same
  subspace and $b$, $Θ(2^(2ℓ))$ pairs per target (32,896 when refuting).
- Candidate arm: a group-law meet-in-the-middle pair table. All $P_i + P_j$ ($i ≤ j$) are
  precomputed once and charged to the factor-base stage; a query walks the $|F| = 256$ base
  points, computes $R − P_3$ and looks it up, for both signs of the target. The table is
  complete, so a miss is a proof.
- Targets: 48 natural (deterministic LCG) + 16 planted-satisfiable (sums of three base points).

== Results (three runs, shared arm64 macOS host, unpinned)

#table(columns: (auto, auto, auto, auto), align: (left, right, right, right),
  [run], [median MITM], [median pairs-and-solve], [ratio],
  [1 (quiet host)], [0.155 ms], [42.764 ms], [276.7x],
  [2 (contended)], [0.124 ms], [11.660 ms], [94.0x],
  [3 (contended)], [0.120 ms], [12.117 ms], [100.7x],
)

- Gate `median ≤ ½ pairs-and-solve`: *PASS* in every run (by $47×$–$138×$).
- Verdict mix: 20 found / 44 refuted of 64; *0 disagreements* with pairs-and-solve on all 64.
- Counted units: at most 1,068 point additions + hash lookups per target (both signs) against
  32,896 quartic solves — sub-$2^(2ℓ)$ in the oracle's own native unit, at $2^(2ℓ)$-entry table
  memory charged to the factor-base stage: 32,896 distinct sums, ~1 MiB retained,
  ~21 ms construction.
- Every witness is re-added in the group and asserted equal to the target (20/20).

#figure(
  image("reports/run-20261015T000000Z/mitm_vs_pairs.png", width: 88%),
  caption: [Per-target wall times of both arms (sorted, log scale) with the half-of-pairs gate line; run `binary_mitm_oracle_l8_20261015`.],
)

== Controls

- `cargo test --release --example binary_mitm_oracle_l8`: at $n = 15$, $ℓ = 5$ the MITM
  decisions agree with the *complete enumeration of all base-point triple sums* (both signs) and
  with pairs-and-solve on 200 random targets, with a found/refuted mix.
- The $ℓ = 8$ run cross-checks every target against pairs-and-solve over the same subspace and
  re-adds every witness in the group.

== Claim boundary

Stage `decomposition` record only. No key recovery, no end-to-end ECDLP, no `vs_rho`, no
asymptotic sub-$ρ$ claim: relation collection remains $Θ(2^n)$ overall because the $2^(2ℓ)$
table replaces the $2^(2ℓ)$ per-target scan. FFD is logged `unmeasured` with reason (group-law
oracle; no polynomial system), which the row's gate permits explicitly. Wall times are
practicality notes from a shared host; the gate metric is the paired median ratio on identical
targets.

== Next rung

$ℓ = 9$ or $ℓ = 10$ on the same design (table $2^(2ℓ)$ entries stays cheap to $ℓ = 10$), a
non-materialized membership-only variant for the factor-base row, and an FFD-logged algebraic
control on the same systems so the SAT/GB frontier comparison carries measured degree data.

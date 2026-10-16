#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.3pt)

= Binary MITM oracle ladder: sub-$2^(2ℓ)$ holds at $ℓ = 9$ and $ℓ = 10$

2026-10-15. `research/binary_mitm_oracle_l9_l10_20261015/`; code `examples/binary_mitm_oracle_l8.rs`
(rung arguments). *Class: engineering* — the $2^ℓ$-fold per-target saving moves with $ℓ$ exactly as
the design predicts; the $Θ(2^n)$ relation-collection total and the $ρ$ frontier are unchanged.

== The row this extends

The binary `decomposition` ledger row landed its first sub-$2^(2ℓ)$ oracle at $ℓ = 8$ (PR #1919,
merge `0f10f333`). That PR's next target was this ladder: the same group-law MITM pair-table design
on fresh frozen instances at $ℓ = 9$ and $ℓ = 10$, holding zero disagreements and the counted
sub-$2^(2ℓ)$ property, with the table still cheap.

== Instances (frozen per rung)

- Curve $y^2 + x y = x^3 + a x^2 + b$, $a = z ∉ F_2$ (non-Koblitz), $b = 1$; deterministic target seed
  fixed in the example; 48 natural + 16 planted targets per rung.
- $ℓ = 9$: $F_(2^27)$, base 519 points over $⟨1, z, …, z^8⟩$; table 134,594 sums, ~4.1 MiB, 601.6 ms.
- $ℓ = 10$: $F_(2^31)$, base 1,053 points over $⟨1, z, …, z^9⟩$; table 554,304 sums, ~16.9 MiB, 4.8 s.
- Reference arm: `SubspaceOracle::low_order(ℓ, 1)` — Gaudry pairs-and-solve on the same subspace.

== Results (shared arm64 macOS host, unpinned, contended by other sessions' builds)

#table(columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto), align: (right, right, right, right, right, right, right, right, right),
  [$ℓ$], [$n$], [verdicts], [disagr.], [max MITM ops], [median MITM], [median pairs], [ratio], [gate],
  [9], [27], [24/40], [0], [2,076], [0.293 ms], [724.9 ms], [2,473.9x], [PASS],
  [10], [31], [20/44], [0], [4,212], [0.735 ms], [3,119.2 ms], [4,244.2x], [PASS],
)

The counted-operation comparison is decisive: 2,076 / 4,212 point additions and lookups against
134,912 / 525,312 quartic solves — the $2^ℓ$-fold gap holds at both rungs and the wall ratio widens
because a quartic solve costs far more than a table lookup. Every witness is re-added in the group
(asserted per target); MITM verdicts agree with pairs-and-solve on all 128 targets.

#figure(
  image("reports/run-20261015T000000Z/mitm_ladder.png", width: 95%),
  caption: [Per-target wall times of both arms at both rungs (sorted, log scale) with the half-of-pairs gate lines; run `binary_mitm_oracle_l9_l10_20261015`.],
)

== Controls

- The ℓ = 8 exhaustive control (complete triple enumeration at $n = 15$, $ℓ = 5$, `cargo test
  --release --example binary_mitm_oracle_l8`) is unchanged and still passes.
- Per-run in-code assertions: every target on-curve, every witness group-verified, every verdict
  cross-checked against pairs-and-solve (0 disagreements, 128/128 targets).

== Claim boundary

Stage `decomposition` records only. The $2^(2ℓ)$-entry table replaces the $2^(2ℓ)$ per-target
scan, so relation collection stays $Θ(2^n)$: no key recovery, no end-to-end ECDLP, no `vs_rho`,
no asymptotic sub-$ρ$ claim. FFD stays `unmeasured` with reason (group-law oracle). Wall times are
L0 practicality notes under external build contention; the gate metric is the paired median ratio.

== Next rung

$ℓ = 11$ (table ~$2^21$ entries, ~68 MiB) and $ℓ = 12$ (the pairs-and-solve wall-clock ceiling of the
reference arm), then the FFD-logged algebraic control on the same systems so the SAT/GB frontier
comparison carries measured degree data.

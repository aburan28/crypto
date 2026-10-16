#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.3pt)

= Binary MITM oracle ladder: sub-$2^(2ℓ)$ holds at $ℓ = 11$ and $ℓ = 12$

2026-10-16. `research/binary_mitm_oracle_l11_l12_20261016/`; code `examples/binary_mitm_oracle_l8.rs`
(rung arguments). *Class: engineering* — the per-target counted saving stays the $2^ℓ$-fold design
factor; the $Θ(2^n)$ relation-collection total and the $ρ$ frontier are unchanged.

== The row this extends

The binary `decomposition` ledger row landed its first sub-$2^(2ℓ)$ oracle at $ℓ = 8$ (PR #1919,
merge `0f10f33`) and held it through $ℓ = 10$ (PR #1944, merge `58da37f`). This round completes the
ladder's next target: $ℓ = 11$ and $ℓ = 12$, where the pair table reaches $2^{2ℓ} = 2^{22}$ entries
and the pairs-and-solve reference arm reaches the ceiling of its historical usability ($ℓ ≈ 12$).

== Instances (frozen per rung)

- Curve $y^2 + x y = x^3 + a x^2 + b$, $a = z ∉ F_2$ (non-Koblitz), $b = 1$; deterministic target seed
  fixed in the example; 48 natural + 16 planted targets per rung.
- $ℓ = 11$: $F_(2^33)$, base 2,017 points; table 2,033,876 sums, 65.1 MB, 2.69 s build.
- $ℓ = 12$: $F_(2^37)$, base 4,023 points; table 8,092,016 sums, 258.9 MB, 14.6 s build.
- Reference arm: `SubspaceOracle::low_order(ℓ, 1)` — Gaudry pairs-and-solve on the same subspace.

== Results (shared arm64 macOS host, unpinned, contended by other sessions' builds)

#table(columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto), align: (right, right, right, right, right, right, right, right, right, right),
  [$ℓ$], [$n$], [verdicts], [disagr.], [max MITM ops], [quartic solves], [median MITM], [median pairs], [ratio], [gate],
  [11], [33], [21/43], [0], [8,068], [2,098,176], [1.746 ms], [1.62 s], [928.0x], [PASS],
  [12], [37], [24/40], [0], [16,092], [8,390,656], [8.566 ms], [18.7 s], [2,187.2x], [PASS],
)

$ℓ = 12$ is the largest rung at which pairs-and-solve was historically usable ("usable wall-clock to
~$ℓ = 12$" was the pre-MITM ledger state); the MITM arm decides the same 64 targets there for at
most 16,092 counted point additions — the $2^ℓ$-fold design gap against 8.4M quartic solves — with
the 259 MB table charged to the factor-base stage. Every witness is re-added in the group (asserted
per target); MITM verdicts agree with pairs-and-solve on all 128 targets.

#figure(
  image("reports/run-20261016T000000Z/mitm_ladder_hi.png", width: 95%),
  caption: [Per-target wall times of both arms at both rungs (sorted, log scale) with the half-of-pairs gate lines; run `binary_mitm_oracle_l11_l12_20261016`.],
)

== Controls

- The exhaustive control ($n = 15$, $ℓ = 5$, complete triple enumeration vs MITM vs pairs-and-solve,
  200 targets, `cargo test --release --example binary_mitm_oracle_l8`) is unchanged and passing.
- Per-run in-code assertions: every target on-curve, every witness group-verified, every verdict
  cross-checked against pairs-and-solve (0 disagreements, 128/128 targets).

== Claim boundary

Stage `decomposition` records only. The $2^(2ℓ)$-entry table replaces the $2^(2ℓ)$ per-target
scan, so relation collection stays $Θ(2^n)$: no key recovery, no end-to-end ECDLP, no `vs_rho`,
no asymptotic sub-$ρ$ claim. FFD stays `unmeasured` with reason (group-law oracle). Wall times are
L0 practicality notes under external build contention; the gate metric is the paired median ratio.

== Next rung

The FFD-logged algebraic control on the same systems (so the SAT/GB frontier comparison carries
measured degree data), and a memory-aware $ℓ = 13$ probe (table $2^{24}$ entries ≈ 1 GiB) if the
factor-base memory gate allows.

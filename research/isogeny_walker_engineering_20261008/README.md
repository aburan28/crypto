# Isogeny-walker engineering preregistrations, 2026-10-08

Two protocols for the repository's isogeny primitives, frozen before any
instrument is built.  Both are **engineering** by §3 of `AGENTS.md`: they
move constants in the walker and the CSIDH toy and touch no security
number, because every prime-field conclusion in the repository rests on
Tate's theorem and the class-level theorems of the structural-completeness
note, which a faster walker does not revisit.  Every measurement is
**PENDING**.

| id | protocol | what it changes | kill line |
|:--|:--|:--|:--|
| K-1 | [precomputed `Φ_ℓ` tables and BMSS kernel polynomials for the P-256 class walker](PROTOCOL-K1-modular-polynomials-bmss.md) | the walker's reach past `ℓ = 61`, and its per-step cost at large `ℓ` | no step below `2×` the current cost at `ℓ = 59`, or any kernel certificate failing |
| K-2 | [√élu in the CSIDH toy for `ℓ` above about 100](PROTOCOL-K2-sqrt-velu-csidh.md) | per-isogeny cost at large prime degree where kernel points are rational | no speedup over plain Vélu at `ℓ = 587` |

Context: the walker (`src/bin/isogeny_walk.rs`, `docs/isogeny-walk/README.md`)
walks the `F_p`-isogeny class of P-256, P-224 and P-192 with
modular-polynomial neighbours and kernel-certified Vélu steps, about 7 ms
per step at `ℓ ≤ 61`, and rebuilds `Φ_ℓ mod p` at every start at a cost
growing like `ℓ⁵` (about 10 s at `ℓ = 59`).  The CSIDH toy
(`src/pqc/csidh.rs`) runs at `p = 419` with plain Vélu and brute-force
point sampling.  Neither is where a security verdict lives; both are
where a drastic constant lives.

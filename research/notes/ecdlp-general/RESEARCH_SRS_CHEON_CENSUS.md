# Cheon's attack on deployed powers-of-tau setups: a census

**Frozen evidence:** `research/srs_auxiliary_census_20260927/` (contract and
sources in its `README.md`, script `srs_census.py`, results in
`results/srs_census.{md,json}`).
**Parent round:** [`RESEARCH_TORSION_AUXILIARY_INPUTS.md`](RESEARCH_TORSION_AUXILIARY_INPUTS.md),
whose §10 item 2 this answers.
**Scoreboard:** `docs/index-calculus-scoreboard.html`, panel
`auxiliary-inputs-20260918`, table *Deployed setups*.

> **Verdict.**  Every deployed powers-of-tau setup examined gives away
> `≈ ½·log₂ q` bits against Pollard rho, where `q` is the largest published
> exponent: from **5.7 bits** for the EIP-4844 mainnet transcript to
> **14.3 bits** for the Perpetual Powers of Tau.  The rule holds because
> `r − 1` on both BLS12-381 and BN254 is smooth enough that a divisor
> within 11% of every `q` exists, so the leak size, not the curve,
> sets the loss.  The published analyses of these setups used the largest
> *power of two* below `q`, which is about `q/2`; the largest admissible
> divisor recovers that factor of two and makes each attack **≈ 0.5 bit
> cheaper than published**.  No setup is broken: the cheapest floor is
> `2^112.3` group operations, and it needs `2^112` memory at that cost.
> Nothing here touches the plain ECDLP.

## 1. What was known, and what this adds

Cheon's `p − 1` algorithm recovers `α` from `G, [α]G, [α^d]G` for any
`d | r − 1` in `√(r/d) + √d` group operations.  A powers-of-tau setup
publishes `[τ^i]G₁` for every `i ≤ q`, so every divisor of `r − 1` up to
`q` is available.  The application to trusted setups is not new: it was
worked through on ethresear.ch (*Cheon's attack and its effect on the
security of big trusted setups*, thread #6692, 2019–2020) for Sapling,
Filecoin, Aztec and the Perpetual Powers of Tau, and Zcash ZIPs issue
#310 tracks it for BLS12-381.

This round reproduces that analysis in this repository's unit and adds
three things:

1. **The EIP-4844 ceremony**, which postdates the thread and is now the
   most widely used setup — every blob on Ethereum mainnet is committed
   under it.
2. **Setup sizes read from each project's own specification**, rather
   than carried over.  They agree with the thread wherever the thread
   states them.
3. **The largest admissible `d`**, not the largest power of two (§4).

## 2. Boundaries

As fixed in the frozen README before computing: reference Pollard rho
with the negation map, `√(πr/4)`; floor `√(r/d) + √d` at constant 1;
unit `log₂` group operations.  **Every cost below is a generic floor, not
a measurement at 255 bits.**  The parent round measured the practical
algorithm (fixed-base tables) at `≈ 20×` its floor over 24–48 bits, which
would add `≈ 4.3` bits; applying that constant here is extrapolation and
is kept out of the headline.

## 3. The census

| setup | curve | largest public exponent `q` | best `d \| r−1`, `d ≤ q` | Cheon floor | rho | **bits lost** |
|:--|:--|--:|--:|--:|--:|--:|
| EIP-4844 KZG, mainnet | BLS12-381 | `2^12 − 1` | `3,648 = 2^6·3·19` | `2^121.51` | `2^127.25` | **5.74** |
| EIP-4844 KZG, largest transcript | BLS12-381 | `2^15 − 1` | `30,531 = 3·10177` | `2^119.98` | `2^127.25` | **7.27** |
| Zcash Sapling | BLS12-381 | `2^22 − 2` | `4,142,391 = 3·11·125527` | `2^116.44` | `2^127.25` | **10.82** |
| Filecoin | BLS12-381 | `2^28 − 2` | `265,113,024 = 2^6·3·11·125527` | `2^113.44` | `2^127.25` | **13.82** |
| Aztec Ignition | BN254 | `100,800,000` | `100,663,296 = 2^25·3` | `2^113.51` | `2^126.62` | **13.12** |
| Perpetual Powers of Tau | BN254 | `2^29 − 2` | `536,259,126 = 2·3·13·29·237073` | `2^112.30` | `2^126.62` | **14.32** |

The `p + 1` variant (needs `2d` powers, `d | r + 1`) is never better on
these setups: `+0.56` to `+5.42` bits dearer than `p − 1`, in agreement
with the parent round's finding that it never wins by more than a few
bits on standard curves.

## 4. The cross-check, and the one correction

| setup | published `d` | this census `d` | security gained by the attacker |
|:--|--:|--:|--:|
| Zcash Sapling | `2^21` | `4,142,391` | `+0.49` bit |
| Filecoin | `2^27` | `265,113,024` | `+0.49` bit |
| Aztec Ignition | `3·2^25` | `3·2^25` | `0` |
| Perpetual Powers of Tau | `2^28` | `536,259,126` | `+0.50` bit |

Every published `d` divides `r − 1` and is `≤ q`, so every published
figure is a valid attack cost; none is wrong.  But the thread's `d` are
the 2-adic divisors — natural, because pairing-friendly curves are chosen
with large 2-adicity for FFTs and the power of two is the conspicuous
factor.  Cheon's algorithm does not need `d` smooth, 2-adic or anything
else; any divisor of `r − 1` serves.  When `q = 2^k − small`, the largest
power of two below `q` is `2^{k−1}`, half of what is available, and the
odd part of `r − 1` supplies a divisor just under `q`.  For Sapling it is
purely odd: `3 · 11 · 125527`.  Aztec is the exception because `3·2^25`
happens to be optimal below 100.8 million.

**Class** (`AGENTS.md` §3): **accounting**.  The algorithm did not change;
the admissible `d` was under-counted.  Half a bit is small, and it is
recorded at its size.

## 5. Reading it

**The leak size is the whole story.**  Exactly, the loss is
`½·log₂ d − 0.17`, the constant being rho's negation-map factor
`√(π/4)`.  And `d/q` ranges `0.89–0.999` over the six setups, so
`½·log₂ d` is within `0.09` bit of `½·log₂ q` on every row.
On these two curves there is nothing curve-specific left to audit: the
number of powers a setup publishes determines what it gives away.  This
is R1 of [`RESEARCH_REPRESENTATION_STRUCTURE.md`](RESEARCH_REPRESENTATION_STRUCTURE.md)
in its plainest form — the same algorithm on the same group is worth
14 bits or 6 bits according to what was published, and nothing else.

**Independent secrets capped the EIP-4844 exposure.**  The ceremony
produced four transcripts, and its participant spec requires *"4
different secrets"*.  The secret behind mainnet blobs is therefore
exposed only through its own `2^12` powers (5.74 bits), not through the
`2^15` transcript (7.27 bits).

**Inheritance.**  A project that took its phase 1 from a larger setup
inherits that setup's `q`, whatever its own circuit size: truncated
parameter files are prefixes of the same `τ`'s powers, and the full
transcript stays public.  Every project built on the Perpetual Powers of
Tau sits on the 14.32-bit row.

**What recovering `τ` would buy differs by proof system.**  For KZG
commitments (EIP-4844, PLONK over Ignition) `τ` is the entire trapdoor,
and knowing it forges openings.  For Groth16 circuits built on a phase-1
powers of tau (Sapling, Filecoin, projects on the PPoT), `τ` is one
component of the trapdoor and the phase-2 secrets are also required.

**Context against the pairing-side attacks — indicative only.**  The
repository's state-of-the-art note puts BN254 at `~100` bits after exTNFS,
so on BN254 Cheon (`2^112`–`2^113.5`) is not the binding constraint.  For
BLS12-381, the thread quotes `117–120` bits for the NFS side; Filecoin's
Cheon floor of `2^113.4`, or `≈ 2^117.8` at the measured constant, sits
at or below that range, Sapling at its edge, EIP-4844 above it.  These
compare a group-operation floor with NFS cost estimates in different
units, and the comparison indicates order only.

**Feasibility.**  The floor's BSGS form stores `√(r/d)` points, `2^112`
at the cheapest row.  The parent round's kangaroo variant removes the
memory for `2–4×` the operations.  Every row is a statement about
security margin, not a break.

## 6. What this is not

- **Not a plain-ECDLP result.**  Every row is handed `[τ^d]G` by the
  setup.  An algorithm seeing only `G` and `[τ]G` still faces `Ω(√r)`;
  the parent round's closing sentence stands.
- **Not a measurement at size.**  The divisors are exact; the costs are
  floors, and the practical constant is extrapolated from 24–48 bits.
- **Not new in kind.**  It extends and corrects a public analysis by
  half a bit and one ceremony.

## 7. Next

- **The G₂ and G_T sides.**  Pairing the G₁ and G₂ powers yields
  `[τ^{i+j}]` in G_T, so `d ≤ q₁ + q₂` is available there, at the price
  of G_T arithmetic.  With `q₂ = 64` for EIP-4844 the gain is negligible,
  but setups with large `q₂` (Sapling, Filecoin, PPoT all publish `2^21`
  to `2^28` G₂ powers) were not priced here.
- **Other auxiliary-input deployments.**  This census covers
  powers-of-tau.  `q`-SDH signature and broadcast schemes that publish
  `[α^i]G` are the rest of the parent round's §10 item 2.

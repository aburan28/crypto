# Cheon exposure of deployed powers-of-tau setups — census, 2026-09-27

Follow-on to `research/torsion_auxiliary_inputs_20260918/`, whose note
(`research/notes/ecdlp-general/RESEARCH_TORSION_AUXILIARY_INPUTS.md` §10,
item 2) names the next step: *protocol audit, not curve audit — which
deployed schemes publish `[α^i]G` for `i` in the thousands or beyond.*
The parent round priced the curves at hypothetical leak sizes (`2^20`,
`2^32`, `2^64`).  This round fixes the leak size to what deployed setups
actually published.

## The question

For each deployed powers-of-tau setup: given the largest exponent `q`
such that `[τ^q]G₁` is public, what is the largest divisor `d` of `r − 1`
with `d ≤ q`, and how far below Pollard rho does that put Cheon's
generic floor?

## Boundaries, stated before computing

| boundary | value | derivation |
|:--|:--|:--|
| reference | Pollard rho with negation map, `√(πr/4)` group operations | the parent round measures `S ≈ 1.3` for it at 24–48 bits |
| floor, DLPwAI, p−1 | `√(r/d) + √d` group operations, constant 1 | Cheon 2006/2010; Boneh–Boyen generic lower bound `Ω(√(p/d))` |
| floor, DLPwAI, p+1 | `√(r/d) + d`, needs `2d` powers, `d \| r+1` | as the parent round's `curve_divisors.py` |
| unit | `log₂` group operations | |

Every cost in the tables is a **generic floor at constant 1, not a
measurement at this size**.  The parent round measured the practical
p−1 algorithm with fixed-base tables at `≈ 20×` its floor over 24–48 bits
(`ops / √(p/d)` = 21.5 / 20.1 / 15.3 / 24.2); applying that constant at
255 bits is extrapolation, and the JSON reports it in a separate field.

## Setup sizes, and where each number comes from

| setup | curve | largest public G₁ exponent `q` | primary source |
|:--|:--|--:|:--|
| EIP-4844 KZG, mainnet transcript | BLS12-381 | `2^12 − 1` | `ethereum/kzg-ceremony-specs` README: four sets, `[τ₁^0]₁ … [τ₁^(2^12−1)]₁`; `docs/participant/participant.md`: *"The participant MUST generate 4 different secrets"* — so the larger transcripts expose a **different** secret |
| EIP-4844 KZG, largest transcript | BLS12-381 | `2^15 − 1` | same; not the transcript mainnet blobs use |
| Zcash Sapling powers of tau | BLS12-381 | `2^22 − 2` | `powersoftau` `src/lib.rs`: `TAU_POWERS_LENGTH = 1 << 21`, `TAU_POWERS_G1_LENGTH = (TAU_POWERS_LENGTH << 1) − 1` |
| Filecoin powers of tau | BLS12-381 | `2^28 − 2` | Filecoin trusted-setup posts: `2^27` τ-powers (64× Zcash); G₁ length `2·2^27 − 1` by the same construction |
| Aztec Ignition | BN254 | `100,800,000` | `AztecProtocol/ignition-verification` `Transcript_spec.md`: `x^1 … x^100,800,000` in G₁ |
| Perpetual Powers of Tau | BN254 | `2^29 − 2` | `privacy-ethereum/perpetualpowersoftau`: *"up to 536870911 powers"* = `2^29 − 1` points |

## Prior work

The same attack on the same setups was analysed on ethresear.ch,
*Cheon's attack and its effect on the security of big trusted setups*
(thread #6692, 2019–2020), with Zcash cryptographers participating;
Zcash ZIPs issue #310 tracks a cost analysis for BLS12-381.  This round
is a reproduction and cross-check, not a discovery.  It adds three
things: the EIP-4844 ceremony, which postdates the thread; exact
setup sizes from each project's own specification; and the largest
admissible `d` rather than the largest power of two (see RESULTS).

## Reproduce

    python3 srs_census.py        # writes results/srs_census.{md,json}

Imports `factor_bounded` and `largest_divisor_below` unchanged from
`../torsion_auxiliary_inputs_20260918/curve_divisors.py`.  Deterministic
(`random.Random(1)`).  BLS12-381's `r − 1` is fully split; BN254's leaves
one 44-digit cofactor unsplit, so its `d` values are exact lower bounds
(a divisor using that cofactor would exceed every `q` here anyway).

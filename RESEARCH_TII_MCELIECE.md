# TII McEliece key-recovery challenges — imported results and ideas

Notes on <https://github.com/mjosaarinen/tii-solved> (Saarinen, pinned commit
`fe91133`), whose public data and verifier are vendored under
[`research/tii_mceliece/`](research/tii_mceliece/README.md). This file records
what the attacks are, what is and is not documented publicly, and where the
ideas touch our own code-based crypto work in `src/pqc/`.

## What the challenges are

TII Track-2 asks for *key recovery* on binary Goppa codes: given only the
public binary parity-check matrix `H`, produce a support
`L = (L_0 … L_{n-1}) ∈ GF(2^m)^n` and a monic irreducible Goppa polynomial `g`
of degree `r` that generate the same code. The 2023 challenge labels ("TII-252"
= claimed `2^252` work factor) are known to be badly inaccurate — the honest
difficulty is driven by `(m, r, n)` and, crucially, by how far `n` is below the
full `2^m`, not by the label.

| Challenge | m | r | n | attack family | wall time | peak memory |
|---|---|---|---|---|---|---|
| TII-254 | 8 | 12 | 223 | certified 64+16 pair-core, GPU | 16 h | 128 GiB |
| TII-252 | 10 | 11 | 1008 | HOVER, `(p,s)=(3,36)` | 6 min 15 s | 2.7 GiB |
| TII-246 | 10 | 11 | 1009 | HOVER, degree-3 | 3 min 11 s | 2.7 GiB |
| TII-240 | 10 | 11 | 1010 | HOVER, degree-3 | 3 min 06 s | 2.7 GiB |
| TII-213 | 9 | 10 | 496 | HOVER, `(p,s)=(4,24)` | 6 h 45 min | ≥144 GiB |
| TII-129 | 9 | 9 | 509 | HOVER, `(p,s)=(4,20)` | 4 h 30 min | 80 GiB |
| TII-248 (control) | 9 | 7 | 482 | HOVER, `(p,s)=(3,17)` | 13 min 38 s | 0.4 GiB |
| TII-83 (control) | 8 | 5 | 253 | MinRank `(p,s)=(4,5)` | 24 s | 1.7 GiB |

Controls TII-83 and TII-248 reproduce keys Hemmert had already published
(ePrint 2026/1339), so they double as a head-to-head: 24 s and 13 min 38 s on a
12-core laptop against 13 min and 17 h 18 min on a 2×128-core server — roughly
30×–75× faster on a small fraction of the cores. That gap is implementation
engineering (Rust/M4RI kernels, tiled outer-parallel elimination), not a
different asymptotic.

## Idea 1 — HOVER: higher-order vanishing distinguishes Goppa from random

Six of the eight recoveries are attributed to *HOVER*, "higher-order vanishing
between binary Goppa and random linear codes". The shape, as far as the
published provenance shows it, is the familiar algebraic-distinguisher line of
attack on Goppa/alternant codes:

- A binary Goppa code with `n ≪ 2^m` is *not* indistinguishable from random —
  certain degree-`p` products/derivatives of codewords (the square-code and
  higher-order analogues) vanish or drop rank in a way a random code of the
  same parameters does not.
- Each recovered key's provenance carries exactly two attack parameters, `p`
  and `s`, plus a seed: `(p,s) = (3,36)` for TII-252, `(4,20)` for TII-129,
  `(4,5)` for the TII-83 control. `p` reads as the vanishing order (the
  degree of the product construction) and `s` as the number of coordinates
  shortened/punctured before it — the standard "shorten to `s`, then look for
  the rank drop" step.
- The run records label TII-129/213 as "degree-`p` MinRank family", tying the
  distinguisher to a MinRank instance: the drop in rank is what a MinRank
  solver exploits to pin down the support, after which `g` follows.
- Cost is then dominated by one exact binary kernel computation. That is why
  TII-129/213 need 80–144 GiB while the `(m,r) = (10,11)` cases finish in
  minutes at 2.7 GiB — memory, not time, is the wall, and the `m = 9`
  instances land on the wrong side of it.

The recoveries are seeded and deterministic (`--seed 25200` reproduces TII-252
end to end), and re-derive **byte-identical** support and Goppa polynomial for
the degree-3 `(10,11)` family.

## Idea 2 — TII-254: pair-core decomposition, not HOVER

TII-254 (`m=8, r=12, n=223`) needed a different algorithm entirely. The
provenance names a "certified 64+16 pair-core decomposition and direct GF(256)
graph-pencil recovery", and the run is structured as:

1. a Krylov sequence on GPU (13 h 36 min of GH200 time),
2. PM-basis (matrix-Padé / block-Wiedemann style basis computation) on 64 CPU
   cores, ~1 hour,
3. eight parallel reconstruction shards, ~58 min each,
4. exact certification, then a Sage finish; once the compact complete-pair
   evidence exists, locator recovery plus full verification takes ~40 seconds.

That is a block-Wiedemann pipeline: Krylov sequence → PM-basis → sharded
reconstruction, with the certificate checked exactly at the end so the GPU
arithmetic never has to be trusted. Total ≈21.3 GH200 GPU-hours, ~16 h critical
wall time, 128 GiB.

The recovered TII-254 key is explicitly an **equivalent** decoding key
(`original_goppa_key_identity_claimed: false`) — it generates the same code but
is not claimed to be the challenge author's own `(L, g)`. For key recovery that
is a complete break; it just means the support/polynomial need not match the
generator's. TII-253 was solved (per upstream) in August 2026; TII-255 remains
open.

The algorithm itself is deferred to a forthcoming publication, so everything
above is read off provenance fields and run records — treat it as a sketch of
the shape, not a specification.

## Idea 3 — the verification discipline is worth stealing

Independent of the attacks, the artifact hygiene here is a good model for our
own `research/` results:

- **Row-space equality, not echelon form.** Accept iff
  `RowSpace(H_rec) == RowSpace(H)` over `GF(2^m)`, where
  `H_rec[j·r+k, l] = (y_l · x_l^k)^(2^j)` with `y_l = 1/g(x_l)`. This is
  independent of how `H` was stored, and of which equivalent key was found.
- **Parameters owned by the verifier.** `(m, r, n)` are checked against a table
  inside the verifier, not read from the candidate, so a key cannot certify
  itself against parameters it chose.
- **Digest binding.** Each key records the SHA-256 of the exact public-key file
  it was recovered against; the verifier recomputes it. `SOURCES.md` pins the
  upstream commits the public keys came from.
- **Canonical JSON, generated views.** The key is stored once as `GF(2^m)`
  polynomial-basis integers plus the field modulus; the TII text and Sage
  pickle forms are deterministic exports. Pickle is flagged as code-executing
  and discouraged — we omit the pickles from our import for that reason.
- **Preconditions gate the expensive check.** Support distinctness, `g` monic
  and irreducible, `g` non-vanishing on the support, and `rank(H) = m·r` are
  all checked before the row-space reconstruction runs.

## Relevance to this repository

`src/pqc/mceliece.rs` and `src/pqc/classic_mceliece.rs` implement the scheme;
nothing here threatens Classic McEliece at its standardized parameters. The
challenges are broken because `n` sits far below `2^m` (e.g. `n = 223` against
`2^8 = 256`, `n ≈ 1008` against `2^10`), which is precisely the regime where
square-code/higher-order distinguishers bite. Standard parameter sets take
`n` at or near `2^m` with much larger `m`, and the published labels being wrong
is a statement about the challenge generator's cost model, not about the
scheme.

Two follow-ups worth doing, neither started:

1. Port the row-space verifier to Rust against our `GF(2^m)` code so the keys
   can be checked without SageMath, and use it as a test vector source for the
   Goppa parity-check construction in `src/pqc/`.
2. Write down the HOVER distinguisher explicitly for small `(m, r, n)` and
   measure at which `n/2^m` ratio the rank drop disappears — that ratio is the
   actual safety margin the standardized parameters are buying.

## References

- `research/tii_mceliece/UPSTREAM_README.md` — upstream results and timings
- `research/tii_mceliece/KEY_FORMATS.md` — field-element mapping, key formats
- `research/tii_mceliece/SOURCES.md` — upstream provenance and digests
- T. Hemmert, ePrint 2026/1339 — the prior solutions the controls reproduce
- `ElenaKirshanova/tii_decoding_challenge` — the official challenge repository

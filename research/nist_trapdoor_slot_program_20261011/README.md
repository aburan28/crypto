# NIST prime-curve trapdoor slot program (2026-10-11)

A preregistered portfolio of experiments, one per place a covert ECDLP
trapdoor could live in the NIST P-192 … P-521 parameters, with one round
executed.

| file | what it is | status |
|---|---|---|
| `PROTOCOL.md` | slot taxonomy, round-0 findings, registered experiments E1a–E5 with predictions, falsifiers and budgets | protocol |
| `census.md`, `census.json` | round 0 output of `examples/nist_slot_census.rs` (35 s, hashes in `PROTOCOL.md` §3) | **executed** |
| `slots.svg` | diagram: slots, invariant vs model-dependent split, status colours | editable source |
| `PROTOCOL.pdf` | `PROTOCOL.md` + `slots.svg` rendered with pandoc/xelatex | rendered copy |

Regenerate round 0:

```text
cargo build --release --example nist_slot_census
./target/release/examples/nist_slot_census          > research/nist_trapdoor_slot_program_20261011/census.md
./target/release/examples/nist_slot_census --json   > research/nist_trapdoor_slot_program_20261011/census.json
```

Rebuild the PDF (pandoc + xelatex; `slots.pdf` is `slots.svg` rendered
with cairosvg):

```text
{ cat PROTOCOL.md; printf '\n\n## Diagram\n\n![Slot taxonomy and round-0 status](slots.pdf)\n\n## Round 0 census (verbatim)\n\n'; sed 's/^# /### /; s/^## /#### /' census.md; } > /tmp/protocol_full.md
pandoc /tmp/protocol_full.md -o PROTOCOL.pdf --pdf-engine=xelatex --resource-path=. -V geometry:margin=2cm -V mainfont="DejaVu Serif" -V monofont="DejaVu Sans Mono" -V fontsize=10pt
```

Round 0 facts worth knowing before reading the protocol:

- all five curves pass every public isogeny-invariant check
  (prime order, cofactor 1, `k > 1000`, no small square factor of
  `4p − t²`);
- the X9.62 procedure fixes `b²`, not `b`; on P-192, P-256, P-384 and
  P-521 the sign of `b` selects between the curve and its quadratic
  twist, and the standard picked the smaller root on three of them and
  the larger on P-384 (P-224 has `p ≡ 1 mod 4`, so no freedom);
- `G` is never the lexicographically first point on any of the five;
- 84 public nothing-up-my-sleeve seed candidates per curve: zero hits;
- the P-224 and P-192 twist orders still have unsplit composites at
  the `2²²` rho budget (pending E5).

## Status update, 2026-10-16: seed provenance resolved

The seeds are now explained (Anthropic Frontier Red Team and F. Valsorda,
2026-10-16; reproduced independently in
`research/nist_seed_provenance_20261016/`, natively in
`src/cryptanalysis/nist_seed_provenance.rs`). Each prime-curve seed is one
20-byte block of `(int)(drand48() * 256)` output after `srand48(t)`, with `t`
a UNIX time between May 1997 and April 1998. `PROTOCOL.md` stays as
registered; this is what the finding does to its slots:

| slot / item | effect |
|---|---|
| S3 seed → `b` | Resolved. The seed is generator output, not a free 160-bit choice. The remaining freedom is the start time and which accepted candidate was published, about 30 bits per curve. |
| S3 root sign | Explained as an algorithm, not a per-curve bit. On every curve with `p ≡ 3 (mod 4)` (P-192, P-256, P-384, P-521 and the X9.62 examples prime192v2/v3 and prime239v1–v3; nine curves, chance 2^−9), `b = (−27/r)^((p+1)/4) mod p`, the root that exponentiation returns. That it is the smaller root on some curves and the larger on P-384 is incidental. P-224's sign carries no freedom (`p ≡ 1 mod 4`). |
| S3 / E3 seed-grinding density | Moot as a test of the NIST curves. The toy-scale densities stay of generic interest. |
| S4 base point | Answered for P-192: `G.x` is the first on-curve `x` from Lenstra's FreeLIP `zrandom` seeded with 20000, after 7 discarded outputs, and `G.y` is the `(p+1)/4` root (E4a can now record this). Open for P-224, P-256, P-384 and P-521: these are not drand48-family output, and for P-384 and P-521 not FreeLIP output under the same recipe for any generator seed. |
| S7 seed as fingerprint | Moot. The 84-candidate dictionary's zero hits are as expected, and no `SHA-1(secret ‖ counter)` family is involved. |
| S1, S2, S5, S6, S8 | Unaffected. |

A new open item comes with the finding: in each 1997–98 run about half of
the candidates that pass every known check were skipped, for a reason no
curve property tested so far explains (`research/nist_seed_provenance_20261016/`).

Nothing here claims or finds a weakness in any NIST curve.  Per
`AGENTS.md`, proposed items stay marked pending until a follow-on PR
lands their evidence.

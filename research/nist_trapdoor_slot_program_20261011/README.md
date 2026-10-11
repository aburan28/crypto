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

Nothing here claims or finds a weakness in any NIST curve.  Per
`AGENTS.md`, proposed items stay marked pending until a follow-on PR
lands their evidence.

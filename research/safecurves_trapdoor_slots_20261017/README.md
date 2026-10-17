# SafeCurves / SECG exotic trapdoor slots (2026-10-17)

Continuation of `research/nist_trapdoor_slot_program_20261011/`: which
trapdoor families could pass the SafeCurves criteria (or the weaker SECG
criteria), and where they would live.

| file | what it is | status |
|---|---|---|
| `PROTOCOL.md` | criteria-versus-trapdoor map (T1-T6), rounds, predictions, follow-ups | protocol |
| `census.md`, `census.json` | round A: sextic-twist-orbit census of secp192k1/224k1/256k1 | **executed** |
| `grind_demo.md`, `grind_demo.json` | round B: toy-scale grind of a strong-twist / weak-sextic-sibling pair | **executed** |
| `twist_orbit_census.py`, `sextic_grind_demo.py` | the two SageMath scripts (SEA point counting, nothing pinned) | code |
| `slots.svg` | diagram: criteria levels vs slots, round-A/B numbers | editable source |
| `rho_costs.svg` | round-A quantitative graph: rho cost per sextic class | editable source |
| `report.typ`, `PROTOCOL.pdf` | rendered report (Typst; protocol + diagram + graph + rounds) | rendered |

Headline facts (measured):

- for `j = 0` curves over `p = 1 mod 6` there are six twist classes with six
  different orders; twist-security criteria certify exactly one (the
  quadratic twist);
- secp192k1 has a sextic class at rho `2^28`, secp224k1 at `2^55.5`,
  secp256k1 at `2^65.5` (the well-known weak secp256k1 twist is a *cubic*
  twist, not the quadratic one the criterion names);
- at 48 and 96 bits, a curve pair passing scaled SafeCurves-style
  published-plus-twist checks with a weak sextic sibling was found within
  2-3 scanned primes - the designer's grind cost is negligible (T5);
- SafeCurves-conforming Montgomery/Edwards models (`j != 0, 1728`) have
  only the quadratic twist, so T5 does not touch them; the slots with no
  public check at any level are T1 (representation) and, in spirit, T3
  (rigidity that is not manipulation-proof).

No deployed-curve weakness beyond the census rows is claimed; no key
recovery. Round C (registered): native Rust rewrite and the 256-bit grind
prediction.

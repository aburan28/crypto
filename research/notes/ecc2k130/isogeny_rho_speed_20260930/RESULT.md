# Result: isogenous rho is not faster than ECC2K-130

**Verdict:** No. Among reachable isogenous curves, Pollard rho is not faster
than on the challenge Koblitz model. Descending leaves lose the cheap
geometric 2-Frobenius, so eligible automorphism order drops from `A = 2n` to
`A = 2`. At `n = 131` that is an `√131 ≈ 11.45×` iteration penalty before any
per-step GPU constant. Modal/RunPod GPU comparison of leaf curves was not
executed: credentials are absent here, and the campaign CUDA stack has no
leaf kernel.

Class of change: **accounting** (AGENTS.md §3) — the headline moves with the
automorphism order that was already the boundary, not with a new algorithm.

## One table, one unit

Unit: `S = charged group operations / √r`. Walk-only diagnostic:
`S_walk = walk_steps / √r`.

### Derived boundary (ECC2K-130)

| row | `A` | `S` floor | ratio to `E_0` |
| --- | ---: | ---: | ---: |
| `E_0` signed-Frobenius (reference) | 262 | 0.07743 | 1.00× |
| descending 263-leaf, negation only | 2 | 0.8862 | **11.45×** |
| plain points (`A = 1`) | 1 | 1.2533 | 16.19× |

Horizontal degree-263 loops return to `j = 1` (`h(−7) = 1`). They do not
supply a second Koblitz model. Source: `evidence/boundary.json`.

### Measured screen (`n = 17`, frozen 273-member census)

| arm | `a₆` | `A` | median `S_walk` | median `S` (total) | verified |
| --- | ---: | ---: | ---: | ---: | --- |
| signed-Frobenius | 1 | 34 | 0.262 | 3.608 | 3/3 |
| negation | 1 | 2 | 1.261 | 1.929 | 3/3 |
| negation | 19167 | 2 | 1.171 | 1.800 | 3/3 |
| negation | 38361 | 2 | 0.566 | 1.097 | 3/3 |
| negation | 57005 | 2 | 0.574 | 1.164 | 3/3 |
| negation | 73353 | 2 | 1.187 | 1.827 | 3/3 |
| negation | 92755 | 2 | 0.961 | 1.613 | 3/3 |
| negation | 107601 | 2 | 1.168 | 1.863 | 3/3 |
| negation | 117287 | 2 | 1.207 | 1.894 | 3/3 |
| negation | 131055 | 2 | 1.093 | 1.769 | 3/3 |

Predicted walk ratio `√(2n) = √34 ≈ 5.83`. Measured median `S_walk`
(non-Koblitz negation) / (signed-Frobenius) ≈ 4.0 on this tiny `r`
(`r = 65587`), same order, high variance. **Total `S` for signed-Frobenius
at this size is setup-dominated** (scalar-mul charge in the table build);
that is an accounting artifact of the toy, not a leaf win. At challenge
size the walk dominates and the `A`-floor controls.

Source: `evidence/n17-rho.json`, `evidence/summary.json`. Frozen census:
`experiments/koblitz_isogeny_cost_sweep.json`.

## GPU (Modal / RunPod)

| check | outcome |
| --- | --- |
| `MODAL_TOKEN_*` in env | unset |
| `RUNPOD_API_KEY` in env | unset |
| `modal` / `runpod` Python packages | missing |
| RunPod GraphQL HEAD | 403 |
| prior repo RunPod receipt | 403 Forbidden |
| leaf GPU rho kernel in `ecc2k130/` | not present (Koblitz-only) |

Receipt: `gpu/access-probe.json`, narrative `gpu/ACCESS.md`. No paid GPU
job was launched. A Koblitz-only Modal smoke would not falsify the leaf
automorphism floor.

## What remains open

- Engineering a generic binary GPU rho for descending leaves (expected
  still `√n` worse in operations).
- Endomorphism walks that use the leaf CM order without cheap coordinate
  squaring — still must beat `A = 262` on `E_0` after charging evaluation.
- IC questions on isogenous models remain with the earlier notes; this
  thread only prices **rho**.

## Reproduce

```sh
python3 research/notes/ecc2k130/isogeny_rho_speed_20260930/derive_boundary.py
python3 research/notes/ecc2k130/isogeny_rho_speed_20260930/gpu/probe_access.py
cargo run --release --locked --example isogeny_rho_compare -- \
  research/notes/ecc2k130/isogeny_rho_speed_20260930/evidence/n17-rho.json
python3 research/notes/ecc2k130/isogeny_rho_speed_20260930/summarize.py
```

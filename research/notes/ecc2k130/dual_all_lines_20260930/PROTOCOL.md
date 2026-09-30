# Protocol: degree-263 dual transport over every rational kernel line

Registered 2026-09-30, before the all-line certificate run. The source is
main at `dc7fa83bc2f6061f6b62927a4adb6a564a83eeae` (which includes
[#991](https://github.com/aburan28/crypto/pull/991)). This extends the eight
saved representatives checked in [#750](https://github.com/aburan28/crypto/pull/750).
It tests an arbitrary-line interface; it does not measure PDP yield or DLP
cost.

## Frozen input and method

Use both archived degree-263 twist-torsion bases, seeds `20260924` and
`20260925`. On the first basis `(U,V)`, enumerate the 264 distinct projective
directions `[1,t]` for `0 <= t < 263` and `[0,1]`, in that order. On the
second independent basis, check the fixed directions `[1,0]`, `[1,1]`,
`[1,2]`, `[1,131]`, `[1,262]`, and `[0,1]`. Derive a
kernel generator `G = U + [t]V` for `[1,t]`, or `G = V` for `[0,1]`. Use `V`
as complement in the first case and `U` in the second. The public Certicom
point `P` is the sign witness; `Q` is the second full-point control. All
coefficients and basis coordinates are canonical. No random stream is used.

The existing `DualTransport` checks exact cyclic order, independent
complement, normalized codomain, and orientation. For **every tested** direction,
verify both `dual(phi(P)) = [263]P` and `phi(dual(phi(P))) = [263]phi(P)` and
the analogous identities for `Q`; verify the infinity cases. Enumerate all
262 nonzero forward and all 262 nonzero reverse *twist* kernel points and
check that they map to infinity. Reject off-curve points even when their
abscissa equals a kernel abscissa. Distinct kernel-abscissa sets must number
264 on the full first basis; all six second-basis fingerprints must belong to
that set.

For each seed, independently reconstruct the forward and reverse generators
and direct full-coordinate Vélu images on `P` and `Q` for the fixed directions
`[1,0]`, `[1,1]`, `[1,2]`, `[1,131]`, `[1,262]`, and `[0,1]` using the
bit-polynomial checker from #703. The two-point reference checks are sampled
deterministically only to bound replay cost; full production dual identities
and all exceptional kernel inputs cover all 270 tested directions, including
all 264 rational kernel lines of the first basis. Include negative
controls for zero/out-of-range/noncanonical line coefficients, invalid basis
points, dependent basis, ambiguous sign witness, and off-curve kernel-
abscissa inputs. Fail closed if a control is accepted.

Frozen input SHA-256:

| Input | SHA-256 |
| --- | --- |
| `research/ecc2k130_direction_review_20260924/twist_torsion_results.json` | `6f1e7b22f3471214d38ec0d3196c88fe9edaf45791e9c05ea67764831fd75368` |
| `research/ecc2k130_dual_transport_20260925/dual_transport.py` | `c913bef64a115e4b42d6675ddbafac10379b166e58ddb94402ebbd44f1e7df9d` |
| `research/ecc2k130_oriented_transport_20260924/oriented_velu.py` | `a1a656cc32efd612250b1ecf50dff18af79a8f62f970ee79841bfd390603423b` |
| `research/ecc2k130_direction_review_20260924/redteam_velu_replay.py` | `c08fbb8de9d54906d783be967a36a6073b2a1d1a3583292badd3adcdd7126790` |
| `research/ecc2k130_relations/fastfield.py` | `b5990fda51700bbfba363251c41f02febfd0fc92edfea4f5ddd647bab472789e` |
| `research/ecc2k130_relations/relations.py` | `0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f` |

## Budget, replay, and decision

Commit and push the interface, certificate, verifier, and a SHA-256 lock for
their exact bytes **before** running this campaign. One local campaign has a
2,100-second wall cap and 2 GiB memory cap. Preserve a `FAIL` or `CENSORED`
receipt with completed directions and the exact failed stage; do not silently
drop a hard line. Separately hash the deterministic certificate fields,
excluding host and timing, and run an independent replay that checks all
identities and exception counts from the frozen inputs. Record setup and
checking time separately; supplied torsion bases are not cold discovery.

`PASS_ALL_LINES` requires 264/264 directions on the first basis and all six
fixed second-basis directions,
all 524 nonzero twist-kernel inputs per direction, every full-point and
reference control, all negative controls, and a byte-stable independent
replay. Otherwise classify the outcome `FAIL` or `CENSORED` and keep the
partial rows. A pass proves only this rational degree-263 dual interface on
the frozen public curve and bases. It does not prove arbitrary extension-field
kernels, endomorphism-ring changes, natural PDP yield, relation rank, or an
ECC2K-130 speedup. Update the canonical scoreboard in the outcome PR.

### Preregistration amendment before source lock or campaign

The first protocol commit `c3d1172a` proposed 264 lines on **both** bases
under 900 seconds. A developmental timing check on `[1,4]` (no campaign
receipt) measured 4.015 seconds for construction and 5.733 seconds including
its controls. Extrapolating 528 lines would exceed that cap by a wide margin.
This amendment retains complete 264-line coverage on one basis and fixed
cross-basis controls, with a 2,100-second cap. It precedes source/input lock
and any all-line result; no observed pass/failure was used to select lines.

# P-256 Dickson-union scalar stabilizer, round 297: result

## Verdict

The exhaustive scalar stabilizer is exactly `{+1,-1}`.  Because negation was
already folded when `FB1hc72514a2a8d3` was built, its effective folded action
has order one: all 164 factor-base logarithm classes remain independent.  The
exact favourable collection boundary therefore stays **13.920747 times rho**.

This is a complete negative result for scalar set-stabilizers of this factor
base, not a small-multiplier sample.  A scalar preserving the 328 signed
nonidentity points must have order dividing both 328 and `n-1`; their gcd is
eight.  The native runner enumerates every one of the eight eighth roots in
`F_n^*`, applies each to every column representative, and finds zero hits for
each of the six roots other than `+1` and `-1`.

The two accepted actions replay on all 328 mappings with zero failures and
form an exact subgroup under four composition checks.  Both independent
invocations emit the same 13,187-byte result.  No relation is reported, the
Round-296 structured-degree gate remains failed, no end-to-end cost is
projected, and no unplanted P-256 relation was attempted.

## Exhaustiveness certificate

The P-256 subgroup has registered prime order

```text
n = 115792089210356248762697446949407573529996955224135760342422259061068512044369.
```

Let `S` be the complete 328-point signed factor base.  If `[a]S=S`, then the
action of `[a]` on `S` is free in cycles of length `ord_n(a)`: for nonzero
`P`, `[a]^j P=P` implies `(a^j-1)P=0`, hence `a^j=1 mod n`.  Therefore
`ord_n(a)` divides `|S|=328`.  It also divides `n-1` by Lagrange's theorem.

```text
divisors(328) = 1, 2, 4, 8, 41, 82, 164, 328
gcd(n-1,328)  = 8
possible orders = 1, 2, 4, 8
```

The least deterministic seed producing an order-eight root is 7.  The runner
enumerates all powers of

```text
r = 0xd9508d5b12903a85e572128007b8bffaf0fa858f0aca47cbe343a4b04b3aee09
```

and verifies eight distinct roots with exact orders covering `{1,2,4,8}`.
This is the entire candidate set licensed by the theorem.

## Root screen

Each row checks all 164 low-y representatives.  A hit means the image's
abscissa is a stored folded column and its ordinate is one of that column's
two exact signs.

| exponent `j` | order of `r^j` | image hits / 164 | unique destinations | stabilizer | accepted replay failures |
|--:|--:|--:|--:|:--:|--:|
| 0 | 1 | 164 | 164 | yes, `+1` | 0 / 164 |
| 1 | 8 | 0 | 0 | no | - |
| 2 | 4 | 0 | 0 | no | - |
| 3 | 8 | 0 | 0 | no | - |
| 4 | 2 | 164 | 164 | yes, `-1` | 0 / 164 |
| 5 | 8 | 0 | 0 | no | - |
| 6 | 4 | 0 | 0 | no | - |
| 7 | 8 | 0 | 0 | no | - |

Every rejected root records the exact first missing affine image and a hash of
the complete 164-image transcript.  The accepted canonical action hashes are

- `+1`: `80bd803fcd7c09ee1dd693dbe96f61e7bd3716849989bcfbde72210f2edc64c5`;
- `-1`: `187bc14858989e9831fb7635ab631d16d634888c5e22a9176031763eb7b89ed9`.

The accepted signed group has order two.  Its induced action on folded columns
has order one and 164 singleton orbits, so direct orbit enumeration gives
`K=164` exactly.

## One boundary, one unit

The unit is favourable P-256 group-addition equivalents divided by `sqrt(n)`.
The rows grant a perfect free decomposition oracle, success probability one,
one operation per sample, and no setup, selector, verification, matrix,
recovery, or representation cost.  They are lower boundaries, not end-to-end
measurements.

| variant | independent log classes `K` | `S` | ratio to rho | correctness | class |
|:--|--:|--:|--:|:--|:--|
| Pollard rho | - | 1.3 | **1.000000** | frozen reference | reference |
| Round-295 unquotiented base | 164 | 18.096972 | 13.920747 | exact favourable boundary | control |
| **Round-297 exhaustive scalar quotient** | **164** | **18.096972** | **13.920747** | exact favourable boundary | **negative** |

The ratio did not move.  This is not an engineering improvement or a local
degree improvement; it closes one proposed source of global logarithm
transport.

## Operations and gates

The complete screen performs 1,312 constant-time scalar multiplications:
335,872 fixed bit rounds, 335,872 group additions, 671,744 group doublings,
1,312 affine normalizations, and 1,312 exact lookups.  Accepted-action replay
adds 328 scalar multiplications, 83,968 additions, 167,936 doublings, and 328
normalizations.  All 164 input columns are distinct, on curve, nonidentity,
and paired with their exact negations.

| gate | status | evidence |
|:--|:--:|:--|
| dependency and input identity | pass | both frozen hashes; 164 canonical sign pairs |
| candidate enumeration exhaustive | pass | free-action theorem; `gcd(n-1,328)=8`; all eight roots |
| every candidate image checked | pass | 8 x 164 = 1,312 images |
| accepted actions replay and compose | pass | 328/328 replays; zero failures; four composition checks |
| scalar quotient reduces boundary | **fail** | folded action order one; `K` remains 164 |
| boundary at or below rho | **fail** | 13.920747 times rho |
| structured residual degree at most 5 | **fail** | Round 296 observes productive degree 6 at 13 variables |
| zero relation FP/FN on complete instances | **fail / unset** | this round runs no relation selector |
| usable relation below `2^103` | **fail / unset** | no complete global solver |
| complete collection below `2^120` | **fail / unset** | no relation-collection projection |
| projected storage below `2^50` | **fail / unset** | no promoted selector |
| promoted | **no** | scalar quotient and prior degree gates fail |

The scalar set-stabilizer route is closed for this exact base.  A successor
must use a non-scalar correspondence, construct a different globally
symmetric factor base without collapsing to a generic rho walk, or solve the
global divisor system by a genuinely lower-degree method.  Repeating bounded
integer-multiplier screens cannot evade this result: every preserving scalar
is already among the eight roots tested here.

## Resources and reproduction

Both native runs were isolated, uncontended, and pinned to reserved CPU 4 on
the AMD EPYC 9V74 host.  The canonical run used 1.082073 s wall, 1.081765 s
user, 0.000000 s system, and 2,372 KiB peak RSS.  The independent run used
1.072995 s wall and reproduced the result byte for byte.

- canonical result: 13,187 bytes, SHA-256
  `53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982`;
- semantic evidence SHA-256:
  `733e9ccc05ae9dafb34aa90932385992b38b2e4531307f1998f165d531f615be`;
- two-run isolation receipt: SHA-256
  `1d29a10bd565fce5c542d28a49b8e8ea56e30a32bb016e6252e0b326bf07be59`;
- independent result replay: byte-identical, SHA-256
  `53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982`.

```bash
cargo test --release --bin p256_union_scalar_stabilizer
cargo clippy --release --bin p256_union_scalar_stabilizer -- -D warnings
cargo build --release --bin p256_union_scalar_stabilizer --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_stabilizer_round297_20261007/isolation.jsonl \
  --label p256-union-scalar-stabilizer-round297-canonical -- \
  target/release/p256_union_scalar_stabilizer \
  --factor-base research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/selected-factor-base.json \
  --round296 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_rr_selector_round296_20261007/rr-selector-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_stabilizer_round297_20261007/scalar-stabilizer-result.json
```

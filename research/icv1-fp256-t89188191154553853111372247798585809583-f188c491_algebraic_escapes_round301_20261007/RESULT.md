# P-256 standard algebraic-escape census, round 301: result

## Verdict

No rho-parity route exists within the registered standard algebraic transfer
families.

The new exact certificate closes the pairing ambiguity.  P-256's subgroup
embedding degree is

```text
ord_n(p)
= 38597363070118749587565815649802524509998985074711920114140753020356170681456
= (n-1)/3,
```

a 255-bit integer.  A direct coefficient representation of `F_(p^k)` has

```text
256*k
= 9880924945950399894416848806349446274559740179126251549220032773211179694452736
bits.
```

The certificate verifies `p^k=1 mod n` and `p^(k/q)!=1 mod n` for every prime
`q` dividing `k`.  This is not a low-degree MOV/Frey--Rück target and is not a
smaller representation of the problem.

The independently recomputed curve invariants also show that P-256 is
ordinary, non-anomalous, has generic `j`, and is defined over a prime field
with no proper smaller finite subfield.  Hash-pinned Rounds 25, 28, 297, and
300 then close their registered boundaries: rational endomorphisms act as
known scalars, the complete registered degree-11 walk supplies no log
quotient, only negation stabilizes the width-164 signed base, and pure cover
sheet multiplicity supplies zero quotient-rank gain.

This is deliberately not a universal impossibility theorem.  A future
P-256-specific construction outside the named families could still exist.
No such construction, structured degree-at-most-five system, executable
recovery, or complete below-rho cost is presently established.  The
width-164 free-oracle boundary therefore remains **13.920747 times rho**; no
unplanted P-256 relation was attempted.

## Requirement-to-evidence status

| requested direction / gate | status | Round-301 evidence or gap |
|:--|:--:|:--|
| pursue fundamentally different factor-base transport | verified within named standard families | exact anomaly, supersingularity, automorphism, pairing, CM, isogeny, scalar, cover, and subfield classifications |
| exact pairing/field-transfer certificate | verified | `ord_n(p)=(n-1)/3`; final and prime-divisor minimality powers checked |
| subgroup preservation and recovery audit | verified where constructions exist | pinned typed results reconciled; absent constructions remain unknown rather than inferred |
| exact deterministic arithmetic | verified | factorization product, trace, discriminant, `j`, order reduction, and imported facts replay |
| applicable non-generic P-256 transfer | not found within boundary | every named family is inapplicable, scalar-preserving, kernel-only, or high-dimensional |
| structured residual degree `<=5` | not established | no surviving P-256 relation construction |
| complete cost at or below rho | failed / unset | optimistic independent-log boundary remains 13.920747 times rho |
| usable relation below `2^103` | not established | no P-256 relation oracle |
| collection below `2^120` | not established | no surviving collector or recovery pipeline |
| storage below `2^50` | not established | no promotable construction |
| actual unplanted P-256 relation | correctly not attempted | promotion prerequisites fail |

The exact Round-301 census is complete.  Achieving parity now requires a new
P-256-specific mechanism outside the named family boundary, not another
parameter adjustment to the screened factor bases.

## Recomputed P-256 invariants

| invariant | exact result | relevance |
|:--|:--|:--|
| field | prime `F_p`, 256-bit `p` | no proper smaller finite subfield |
| group order | prime `n`, cofactor 1 | one prime cyclic rational group |
| trace | `89188191154553853111372247798585809583` | `n != p`; not anomalous |
| Frobenius discriminant | `-455213823400003756884736869668539463648899917731097708475249543966132856781915` | exactly matches Round 25 |
| ordinary / supersingular | ordinary / no | prime-field trace is nonzero |
| `j` | `7958909377132088453074743217357398615041065282494610304372115906626967530147` | neither 0 nor 1728; no exceptional automorphism case |

The `j` value is also the value embedded in the ICV1 identity record; this
round recomputes it from `a` and `b` rather than trusting the slug.

## Exact embedding-degree certificate

Round 25's certified factorization is re-multiplied to `n-1`.  Repeated
prime-order reduction succeeds only once, for the factor 3.  It fails for 2
and every other prime factor.  Thus

```text
k=(n-1)/3
 = 2^4 * 71 * 131 * 373 * 3407 * 17449 * 38189
   * 187019741 * 622491383 * 1002328039319
   * 2624747550333869278416773953.
```

The final certificate contains one non-one residue for `k/q` for each of
these eleven distinct prime divisors.  It performs 25 counted modular
exponentiations:

| counted operation | exact count |
|:--|--:|
| modular exponentiations | 25 |
| modular multiplications | 2,881 |
| modular squarings | 5,873 |
| factor divisibility tests | 12 |
| exact integer divisions | 23 |
| integer multiplications | 17 |

These counts measure certificate generation, not an extension-field DLP.
No finite-field runtime or asymptotic crossover is extrapolated from them.

## Named-family census

| family | exact applicability result | quotient / destination result | status |
|:--|:--|:--|:--|
| Smart anomalous curve | `n != p` | required map does not apply | closed |
| supersingular shortcut | curve ordinary, trace nonzero | required hypothesis fails | closed |
| exceptional automorphisms | `j != 0,1728` | only generic negation | closed |
| MOV/Frey--Rück | exact `k=(n-1)/3`, 255 bits | enormous extension representation, not low-degree | closed as low-degree route |
| CM/GLV endomorphisms | minimum noninteger degree about `2^255.975`; rational action scalar | no non-generic rational action | Round 25 boundary closed |
| degree-11 isogenies | 4,096 certified edges / 4,097 models | coordinate transport, no log quotient | Round 28 registered walk closed |
| width-164 scalar stabilizer | all eight possible eighth roots screened | folded action order 1; `K=164` | Round 297 closed |
| cover/Jacobian sheets | 374,400 complete lifted rows | maximum quotient-rank gain 0 | Round 300 pure-fibre family closed |
| base-field subfield descent | `F_p` is prime | no proper smaller finite subfield | closed |

The isogeny statement is scoped to Round 28's complete **registered walk**, not
the entire mathematical isogeny class.  The cover statement is scoped to
kernel-sheet multiplicity; a construction that makes distinct collapsed
decompositions cheaper remains an explicit open obligation.

## Boundary table

The unit is favourable P-256 group-addition equivalents divided by `sqrt(n)`.
Candidate rows omit solver, verification, collection, sparse linear algebra,
and recovery costs and are lower boundaries, not measured attacks.

| variant | optimistic / rho | complete cost measured here | status |
|:--|--:|:--:|:--|
| Pollard rho | **1.000000** | reference | matched boundary |
| minimum-width `FB1hc72514a2a8d3` | **13.920747** | no | unchanged after named-family census |
| registered 17-term `FB1h2f8621cda105` | 394.425280 | no | unchanged comparison |
| future construction outside census | unset | no | open; no algorithm or projection |

This round supplies no new numerator that could be compared with rho.  It
classifies the named transfer denominators and retains the existing negative
boundary.

## Typed transfer assessment

The companion assessment types the source, destination, map, field,
dimensions, subgroup, kernel/recovery obligations, evidence scope, and weakest
unresolved obligation for each family.  Its narrowest supported conclusion is:

> No applicable rho-parity route exists within the registered anomalous,
> supersingular/pairing, exceptional-automorphism, CM/GLV, degree-11 isogeny,
> scalar-stabilizer, pure cover-fibre, or base-field subfield families.

The weakest open obligation is a new construction outside this list that
produces distinct collapsed equations or an easier destination, with
executable inverse recovery and complete one-target cost below rho.

The transfer profile's main workflow was available.  Its referenced
methodology and JSON-template resources were unavailable, and the assessment
records that limitation.

## Gates

| gate | status | evidence |
|:--|:--:|:--|
| dependency hashes and schemas | pass | Rounds 25, 28, 297, 300 pinned exactly |
| curve invariants recomputed | pass | P-256 `p,a,b,n,h`, trace, discriminant, `j` |
| `n-1` factorization | pass | exact product |
| embedding degree | pass | final power plus all prime-divisor minimality checks |
| imported certificates reconciled | pass | four source boundaries match exact fields |
| applicable non-generic transfer found | **fail** | none within named families |
| explicit new recovery | **fail / unset** | no surviving construction |
| structured residual degree at most 5 | **fail / unset** | no surviving system |
| parity at or below rho | **fail** | lower boundary still 13.920747 times rho |
| per-row, collection, storage gates | **fail / unset** | no candidate pipeline |
| promoted | **no** | parity-critical gates fail |

## Resources and reproduction

The canonical and independent runs were pinned to reserved CPU 4 on the AMD
EPYC 9V74 host and were uncontended.  The canonical run used 0.010914 s wall,
0.007154 s user, 0.003651 s system, and 5,960 KiB peak RSS.  The independent
run used 0.018336 s wall and 5,992 KiB peak RSS.  The result and assessment
artifacts are byte-identical.  Wall time is secondary; exact arithmetic and
hash certificates are the evidence.

- result: 16,618 bytes, SHA-256
  `ae4b44315181c18a40a790ae1bc784ffdb05e3061d442b8f1be78a1f25dc20e9`;
- result semantic evidence SHA-256:
  `69982b49633ddb1ca8226360a4797a4d6c80562ae36310c6a6e9545e3dd0c141`;
- transfer assessment: 10,088 bytes, SHA-256
  `4d9168a57f0c848abecf37ff721edd9be30ccb444948f42de0343384b94caeea`;
- assessment semantic evidence SHA-256:
  `709fdc4d0c68f60651454a705f59aefddc82a1380271d10b5216b65bd69e2cd6`;
- two-run isolation ledger: 4,949 bytes, SHA-256
  `6421b20035963d2b6ef4912f9a0addc0443c98f9eb525c317b416d428ce58f94`.

```bash
cargo test --release --bin p256_algebraic_escape_census
cargo clippy --release --bin p256_algebraic_escape_census -- -D warnings
cargo build --release --bin p256_algebraic_escape_census --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label p256-algebraic-escapes-round301-canonical-v1 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_escapes_round301_20261007/isolation.jsonl -- \
  target/release/p256_algebraic_escape_census \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round28 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_isogenous_dickson_round28_20261006/isogenous-dickson-result.json \
  --round297 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_stabilizer_round297_20261007/scalar-stabilizer-result.json \
  --round300 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_cover_fiber_round300_20261007/cover-fiber-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_escapes_round301_20261007/algebraic-escape-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_escapes_round301_20261007/transfer-assessment.json
```

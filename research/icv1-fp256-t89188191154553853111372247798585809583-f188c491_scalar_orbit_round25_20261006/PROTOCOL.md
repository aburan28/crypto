# P-256 exact-endomorphism and scalar-orbit factor-base screen, round 25: protocol

Date frozen: 2026-10-06

Round 24 isolated the only remaining collision-level parity target: exact
elliptic-log transport reducing a geometric factor base to one log class,
combined with a global low-cost selector that is not Pollard rho under another
name.  This round screens the complete P-256 endomorphism ring and the finite
scalar orbits permitted by the prime subgroup order.  It may construct an exact
one-class factor base; it may not call that an index-calculus improvement unless
membership, selector, anchor, and end-to-end costs also pass.

This is a bounded factor-base and obstruction screen, not a P-256 discrete-log
attack.

## Frozen curve and dependency

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- relation model: signed, distinct-column `S17`, balanced `8+9` split;
- comparison base: `FB1h2f8621cda105`, 131,458 columns;
- useful-row comparison target: 138,031;
- round-24 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/negation-result.json`;
- required dependency SHA-256:
  `d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c`.

## Exact CM-order screen

Recompute `t=p+1-n` and `Delta=t^2-4p`.  The following candidate
factorisation is discovery input only:

```text
|Delta| = 3 * 5 * 456597257999
          * 1428624589419343516204097
          * 46523541035814968339936406074986559003387.
```

Do not trust the source of the factors.  Verify their exact product and prove
each factor prime.  Values fitting `u64` use deterministic Miller--Rabin;
larger values require replayable Lucas/Pocklington certificates with a fully
factored `q-1`, including a recursive certificate for
`1869236796843064056413`.  Record every witness and factorisation used by the
certificate.

If all factors are distinct primes and `- |Delta|` is a fundamental
discriminant, conclude exactly that the Frobenius order `Z[pi]` is maximal and
therefore equals `End(E)`.  For `D=Delta congruent 1 (mod 4)`, independently
verify the norm formula

```text
N((a+b*sqrt(D))/2) = (a^2 + |D|*b^2)/4,  a congruent b (mod 2).
```

Report the exact minimum degree `(1+|D|)/4` of a non-integer endomorphism.
Also reduce every `a+b*pi` action on `E(F_p)` using `pi(P)=P`; any resulting
known scalar action is classified as generic scalar transport, regardless of
the algebraic degree of the endomorphism.

## Exhaustive scalar-orbit-width screen

Verify the complete factorisation

```text
n-1 = 2^4 * 3 * 71 * 131 * 373 * 3407 * 17449 * 38189
      * 187019741 * 622491383 * 1002328039319
      * 2624747550333869278416773953.
```

Prove the large final factor prime with a replayable certificate, and prove
`n` prime by Lucas using the complete factorisation.  Choose the least integer
primitive root modulo `n`, recording and replaying every exact-order check.

Enumerate every divisor `r | n-1`.  A scalar of odd order `r` produces `r`
negation-folded columns; a scalar of even order produces `r/2` columns because
the unique order-two element is `-1`.  Rank every attainable width by absolute
distance from 131,458, with smaller width then smaller order as deterministic
ties.  Select the closest width, derive an exact-order scalar, and exhaust all
generators of that cyclic subgroup to choose the representative with minimum
signed binary double-and-add cost.  Report the complete divisor count and the
top width candidates; sampled orders receive no exhaustive label.

## Native factor base and replay

Derive an anchor by SHA-256-labelled try-and-increment lift to P-256, taking
the lower square root and recording no construction scalar relative to `G`.
Construct the selected orbit `P_i=[m^i]H`, fold negation, sort columns by the
repository 33-byte point key, and form the canonical
`ecbench.factor_base/v1` FB1 preimage.  The optional full dump must use
`ecbench.factor_base_dump/v1-wide`, materialise both signs, and reproduce the
same FB1 identity and point digest on rebuild.

Replay every claimed adjacent transport relation in the P-256 group and
verify the wrap/negation edge exactly.  Independently check coefficient order,
all point keys, curve membership, nonidentity, uniqueness up to sign, and the
wide-dump rows.  Record group operations, logical RAM, dump bytes, replay
failures, false positives, and false negatives.

## Information and end-to-end accounting

Separate three cases:

1. **unknown anchor**: all base logs are known multiples of one unknown `h`;
   ordinary zero relations reduce to coefficient identities and provide no
   information about `h`; at least two independent target equations are needed
   to eliminate `h`;
2. **known anchor `H=G`**: all base logs are already known, so the remaining
   target representation is a generic collision/claw search and must be
   labelled rho-equivalent;
3. **coordinate/summation-polynomial use**: record the exact support-polynomial
   degree and multiplication-map degree.  A short scalar evaluation chain or
   a one-class log quotient is not a measured low degree of regularity.

Use round 24's exact negation-folded cross-colour law, including the exact
distinct-column probability:

```text
E[T_K] = sqrt(2*n/p_disjoint(B)) * Gamma(K+1/2)/Gamma(K).
```

Price `K=1` for a known anchor and `K=2` for an unknown anchor.  Report the
optimistic oracle cost, direct eight-addition sample cost, one-list storage,
and ratios to `1.3*sqrt(n)`.  No ordinary relation-row savings may be counted
when all rows are predetermined scalar identities.  Sparse linear algebra is
zero only in the already-known-log case, not evidence of an IC speedup.

## Promotion and stop conditions

Promotion requires all of:

- zero false positives and false negatives on complete checked instances;
- exact replay of every reported transport relation;
- measured structured residual degree no greater than 5;
- projected relation collection below `2^120` operations and `2^103` per
  usable row;
- peak projected materialised storage below `2^50` bytes;
- complete one-target cost below matched rho; and
- a selector and anchor treatment demonstrably outside generic rho/claw
  search.

Attempt no full-depth unplanted P-256 relation unless every gate passes.  If
the orbit reaches one log class but fails membership, anchor, selector, or
end-to-end gates, publish that split result and the exact obstruction.  Emit
deterministic canonical JSON with hashes, a result report, dashboard update,
tests, and an optional native wide dump.

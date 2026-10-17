# Primary-source review and comparator contract

No new external timings were run or quoted as reproduced.

1. [Current native group-law rootfinding source](https://github.com/aburan28/crypto/blob/main/src/cryptanalysis/ic_framework/rootfind.rs),
   blob af158961296bb48ed400e5779825c1710c69d424, introductory contract:
   rootfind-mitm stops at the first hit and acts by the group law, never by
   an algebraic system. It is a strong rational-decomposition competitor
   within that contract. It is not an all-original-polynomial-roots baseline.
   Lookup work is separately identified in that source; preserve its boundary.
2. [Generic wide-prime F4 study](https://github.com/aburan28/crypto/blob/main/research/f4_generic_field_20261015/README.md),
   blob a3c14ae47d7ff670a2d342778b083e34e5fae8fd: genuine generic F4 backends
   and P192/P224/P256 tower workloads exist. Do not repeat a blanket claim
   that wide-prime F4 arithmetic is absent. That study calls its costs stage
   diagnostics on a shared host; basis/decomposition checks do not by
   themselves establish our complete all-original-root output contract.
   Export, auxiliary fibers, representation and complete extraction still
   require a matched adapter. P384/P521 coverage must be checked, not inferred.
3. Cache C1 is a small fixed-degree Frobenius operator, a classical semilinear
   construction. This review establishes no priority or novel complexity.
4. [Identity certificate format](https://github.com/aburan28/crypto/blob/main/docs/identity-certificates/README.md)
   fixes modulus2^61-1. It cannot directly certify an identity relying on
   2=0: (a+b)^2=a^2+b^2 is false over that odd field. A valid route is
   to certify the integer identity with its 2*cross-terms, then formally
   discharge the characteristic-two and quotient-ring implications.
   No certificate ID or Lean-checked claim is fabricated here.
5. The next polynomial/field comparison must include the fastest available
   FLINT/NTL/Magma implementation of the same computation and basis.
   The final solver gate retains genuine matched all-roots F4/F5/F6,
   including strongest available external F4/F5 configurations. Absence
   of an external tool leaves comparison incomplete; it is not a win.

This is a scope and correctness review of identified primary code and its
authors' report, not an exhaustive literature novelty audit.

# Wide odd-field supplied-quintic validation

The complete mathematical contract and proof of the frozen v1 root adapter
remain. This version adds an all-field-root certificate to the filtered
interval output: every distinct field root/multiplicity and the monic
G=gcd(H,X^q-X). Multiplying the distinct output linear factors must equal G.
Every field root is checked in H and exact division recovers multiplicity.
Interval roots are exactly those canonical coordinates below2^d. A capped
split returns unresolved, false component-completion and null exact counts.

No randomized root guesses are accepted without exact coverage. The
factorization step is classical odd-characteristic equal-degree factoring;
no novelty, cryptographic security, generalized-parameter-search or8x claim.
Generic q rather than only characteristic p makes the proof apply to odd
extensions; concrete backends here are prime fields and odd quadratics.

Freeze66 supplied monic quintics across11 fields before the wide run:
Fp3/5/7/11; P192/P224/P256/P384/P521; GF9/GF25. Six patterns include
inside/outside interval roots, double and quintuple multiplicities, two
double factors and an irreducible quadratic squared times a linear factor.
The irreducible quadratic X^2-c is certified by c^((q-1)/2)=-1; no field
root is silently discarded and repeated nonrational factors remain absent
from field-root output with their narrower algebraic meaning explicit.
Expected results come from literal root bags or the irreducibility witness,
not from the splitting algorithm. These are known-answer component tests,
not original hard-system benchmarks or independent external factorization.

3624 small controls rerun the independent Hasse multiplicity oracle before
freeze. Wide-field input seed2026101761, interval widths40/47/53/79/106 for
the five NIST primes, source-identified preset characteristic constants.
No prime interval enumeration or false proper-subfield assumption is used.
90seconds/4GiB/oneCPU. Full serialization is checked before acceptance and
late counts/flags are demoted. The supervisor charges the full cold lifetime.
Output includes diagnostics only; no completed matched F4/F5/F6 or fastest
external comparison, A/A panel or qualified isolation is asserted.

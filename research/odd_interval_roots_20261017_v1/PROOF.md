# Exact odd-characteristic supplied-polynomial interval roots

This is a classical factoring adapter, not a new generalized solver or a
novelty claim. Mathematical contract: all distinct roots in K of a supplied
monic H, filtered by canonical-coordinate interval [0,2^d), with exact
multiplicities. K may be an odd prime or extension field; the current concrete
backends are the existing prime and odd-quadratic implementations. General
extension-field adapters and original multivariate search are not implemented.
Only small correctness controls are prepared in this version; no wide-field
experiment or performance comparison is executed or claimed.

Let q=|K|. X^q-X has exactly the elements of K as simple roots, since its
derivative is-1. Compute G=gcd(H,X^q-X) using modular exponentiation before
Euclidean GCD. G is exactly the monic squarefree product of every field-root
factor, even if H is inseparable or has repeated nonrational factors.
For squarefree split G, standard odd-characteristic equal-degree splitting
uses gcd(G,A^((q-1)/2)-1). Random polynomial A separates field roots by their
quadratic character. Every successful proper split is exact, including the
roots at which A vanishes. Recursively divide into linear factors.
A capped64-attempt branch returns unresolved rather than an empty root set.
A future worker must demote exact counts/completion on this unresolved outcome.

Before accepting output, multiply the distinct linear factors and require
the result equal G. This certifies root coverage given the exact GCD step;
evaluate every root in H. Recover multiplicities by exact linear division
of the ORIGINAL H, then filter interval membership. Filtering never scans
the interval or assumes a proper subfield subspace in Fp. Nonrational factors
are correctly absent from K-root output; they are not declared rational.
Interval roots and all-algebraic-roots-in-domain support are distinct contracts.

Independent controls exhaust all monic quintics over F3/F5 and sample GF9/GF25,
using every tiny field point and Hasse derivative multiplicities.3624 cases
pass, including repeated/inseparable and nonsplit inputs. No sqrt shortcut or
binary-specific characteristic formula is transferred to odd fields.
The complete prime wide-field panel P-192/P-224/P-256/P-384/P-521 with fresh
inputs, exact transcript/output verification, cold costs and bounded failures
is the next registered experiment, not a completed result of this file.

# Frozen degree-three equivariant-coordinate identity test

2026-10-13 local. This test follows the completed [orbit support screen](REPORT.md)
and precedes execution of the new coordinate verifier. It uses the same
four `(p,r)` pairs, `E:y²=x³+1`, `T=(0,1)`, quotient `phi`, and destination
curve as [PROTOCOL.md](PROTOCOL.md).

## Hypothesis

Define the projective coordinate `chi=(y+1)/(y-1)` on `E`, with
`chi(O)=1`, `chi(T)=infinity`, and `chi(2T)=0`. Define the Möbius action
`f(z)=1/(1-z)` on `P¹`, with `f(infinity)=0` and `f(1)=infinity`.
The proposed identity is `chi(P+T)=f(chi(P))` for every rational point,
including its three exceptional kernel points, and `f³=id`. Since the
three points in `ker(phi)` have three different `chi` values, `chi` does
not factor through the quotient isogeny. This is a concrete candidate
for a torsion-equivariant decomposition coordinate, distinct from
the invariant coordinate `x∘phi`.

## Test and limits

Use native Rust arithmetic independent of the prior measured runner.
Enumerate all `p+1` rational source points on each of the four fields.
Check source addition, every curve equation, `T` order three, the
coordinate identity for `P+T` and `P+2T`, the projective action order for
all `p+1` points of `P¹(F_p)`, and the three distinct `chi` values on the
kernel. Check `phi(P+T)=phi(P)` as the invariant control. Abort and retain
any counterexample. Record exact commands, source/binary hashes and raw
counts. The algebraic derivation and proof map will be in the report.

This identity test measures no decomposition, relation solve, or runtime
gain. Passing it licenses a *candidate coordinate* for the pending
extension-field oracle test only. The published comparison is
[Faugère et al., EUROCRYPT 2014](https://www.iacr.org/archive/eurocrypt2014/84410158/84410158.pdf):
it distinguishes nontrivial equivariance from coordinates that factor
through an isogeny and requires the polynomial and complete-DLP costs
to be charged separately.

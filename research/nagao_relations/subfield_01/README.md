# Binary subfield curves and larger factor bases

This successor extends the hybrid to
\(E: y^2+xy=x^3+A x^2+B\), with \(B\ne0\), including coefficients in
\(\mathbb F_4\setminus\mathbb F_2\). It uses even-degree ambient fields
\(\mathbb F_{2^{18}}\) and \(\mathbb F_{2^{30}}\). This is a functionality
extension and a solver-stage engineering diagnostic. It does not claim a
full ECDLP speedup, an exponent improvement, or the original 20% cost goal.
The earlier curve was already defined over F2; merely renaming that case
as a subfield curve would not have extended its scope.

The frozen contract increases the dimension of V from 6 to 7 to 8:
64, 128, and 256 possible abscissas before removing zero and testing curve
membership. V remains an **F2-linear** ONB coordinate subspace. It is not
the coefficient field, an F4-linear base, or the subgroup E(F4).

## Exact extension

For R=(r,s), take
\[
 f=X^2+aX+c+bY,\qquad c=r^2+ar+b(r+s),\quad b\ne0.
\]
The involution sends Y to Y+X, so the norm is
\[
 N=X^4+(b+b^2)X^3+(a^2+ab+A b^2)X^2+bcX+c^2+B b^2.
\]
The pinned value makes N(r)=0. Writing N=(X+r)H gives
\[
 h_2=b^2+b+r,\quad h_1=a^2+ab+A b^2+r h_2,\quad h_0=bc+r h_1.
\]
Every allowed decomposition supplies this normalized function: its zero
divisor consists of the three factor points and -R, with poles 4O. The
Riemann--Roch basis is {1,X,Y,X²}. A function without the X² coefficient
has pole degree at most three, so cannot have those four distinct zeros.
The b=0 case has a squared quadratic norm and cannot have four distinct
abscissas. Thus the normalization loses no decomposition in the declared
distinct-abscissa chart.

The support argument is unchanged. Enumerate h2 in V and solve
b²+b=h2+r. For z in V excluding 0,r,h2, the equation H(z)=0 is
\[
 (r+z)a^2+bz a+K=0,\qquad K=H(z)|_{a=0}.
\]
Set t=bz/(r+z) and a=tw. This gives
w²+w=K/((r+z)t²). Both quadratic steps are solved through a precomputed
linear image of w↦w²+w on the full ambient field, with kernel {0,1}.
Unlike the predecessor's half-trace implementation, this works in even
degree too. No scan over all ambient-field coefficients is introduced.

Once a candidate has H(z)=0, put u=h2+z and v=h1+zu. Then
H=(X+z)(X²+uX+v). The residual quadratic has its distinct roots in V
exactly when v is in the image of w↦w²+uw restricted to V. Its kernel is
{0,u}, so its rank is d-1. Gaussian elimination stores preimages; a
successful membership test recovers w and w+u. The exclusions H(0)≠0,
H(r)≠0 and distinctness of all three roots make this equivalent to
H dividing the squarefree subspace polynomial L_V in the restricted chart.
Extract Y=(X²+aX+c)/b, check the curve equation, and verify the signed sum.

The implementation tests pivot bits directly in the canonical internal
field representation, avoiding repeated full ONB coordinate conversion
inside Gaussian reduction. Every image/preimage update is still a counted
field addition. Image setup and batch inversions are charged cold per cell.
The conditioned branch count remains O(2^(2d)); this is a hybrid with one
chosen abscissa, not the original fully root-free solver circuit.

## General Semaev controls

Both controls use the same non-F2 B, bases, targets, domain, extraction,
and exact signed verification. A does not appear in the abscissa polynomial;
it does affect which abscissas lift on the selected curve.
\[
 S_3(x,y,z)=(xy+xz+yz)^2+xyz+B.
\]
For elementary symmetric e1,e3,e4 in x1,x2,x3,r,
\[
 S_4=e_3^4+(e_4+B)e_3^2+e_4(e_4+B)e_1^2+B^2e_1^4.
\]
`semaev.validateIdentity()` verifies this identity exactly in
F2[x1,x2,x3,r,B] against the quadratic resultant of two S3 polynomials.
The identity is independently checked algebraically here; the B=1 frozen
control is not silently reused. Both formulations compile to the existing
ONB Boolean/XOR circuits and run CryptoMiniSat with one thread. Algebraic
models with no signed point lift are rejected and their projections blocked.
The transformed B=1 S4 chart is not claimed as a general-B control.

## Boundary, falsification, and accounting

Before measuring: with M signed factor points and three distinct abscissas,
there are at most 8 binom(M/2,3) signed triples. Uniform affine target support
is at most min(1,8 binom(M/2,3)/(#E-1)), before the target exclusions. Growing
d moves this counting boundary; it must not be reported as a gain against
a fixed boundary. The matched references are chained S3 and symmetric S4.

Any exact-set or certificate mismatch falsifies correctness. All eight
dimension-eight targets must enumerate within the supplemental 60 seconds
to establish feasibility on this panel. The common five-second trials
remain a separate comparison; supplemental time is never pooled into them.

The repository's parent accounting contract and AGENTS benchmark gate apply.
The frozen WDSat Koblitz corpus does not encode these curve coefficients or
even-degree fields. This is an explicitly exploratory matched stage panel,
not an equivalent completion of the required 60-input, three-repetition
regression, and not a promotion-gate run. No full-DLP or calibrated speedup
is claimed. Full-DLP S, rho and floor ratios, common operations, and speedup
remain null. The field API counts are useful only within the hybrid; the
SAT controls lack comparable complete operation counts.

Targets include a separate uniform affine stratum and a diagnostic stratum
generated from signed triples at d6. Development and holdout seeds are frozen
before execution. Each target is retained at all three dimensions and for
every solver. First verified relation and complete enumeration have separate
records. Setup, search, extraction and verification use exclusive timers;
all failed attempts remain in the totals. Harness target generation and
independent oracle costs are recorded separately and are not supplied to
the solvers. No relation-matrix independence is asserted.

## Reproduction

From the repository root, with pycryptosat installed:

```sh
python research/nagao_relations/subfield_01/run.py --validate-only
python research/nagao_relations/subfield_01/run.py
```

The campaign refuses to overwrite raw evidence. Copy the experiment to a
new sibling directory for a rerun, keeping its imports at the same depth.
Its raw records contain source/configuration hashes, exact expected sets,
point certificates, cold phase costs, and explicit timeouts. Results and
the canonical scoreboard are appended after the frozen campaign finishes.

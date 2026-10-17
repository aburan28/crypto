# Exact non-Q-root repeated tangency branch

For the actual binary curve and target, D has dimension18, offset0, and
z lies outside D. Work on the same four fixed h,k norm families; this does
not enumerate all h,k or solve squarefree norm candidates.

Let x be a repeated source abscissa with x!=0,z and Q(x)=hx+k!=0.
A norm root forces y=P(x)/Q(x) in K to satisfy the curve; hence nonrational
curve fibers yield no candidate here. Because x!=0 the two curve branches
are unramified and distinct. Norm multiplicity>=2 on the vanishing branch
is equivalent to the tangent derivative of F=P+QY vanishing. Implicit
curve differentiation gives Y'=(Y+X^2)/X and therefore
x*F'=(h+1)*x^3+k*x^2+x*s+k*y.
Writing w=y/x gives exactly
s=(h+1)*x^2+k*(x+w),
r=[s*(x+z)+x^3+z^3+(h*x+k)*y+(h*z+k)*yT]/(x+z)^2.
For every rational curve point and admitted x these are necessary and
sufficient and uniquely determine r,s. Independent exhaustive tiny checks
compare the resulting sets to literal H(x)=H'(x)=0 over all r,s.

The curve ordinate equation becomes w^2+w=x+a+b/x^2. A binary Gaussian
basis of the Artin-Schreier map supplies either no solution or one w;
then w+1 supplies the other, exactly. No half-trace parity restriction is
assumed. The two-point AS kernel is precisely {0,1} in a field.
Controls compare all28 small-field right-hand sides to literal enumeration.

At x=0,Q(0)=k!=0, the norm's constant coefficient is P(0)^2+k^2*b and
its linear coefficient is k*P(0). As b!=0, they cannot both vanish.
Thus repeated x=0 belongs only to the Q-root branch handled separately.
The worker excludes the pole because z is outside D; general pole domains
need an additional branch. h=0 cases are not in these fixed profiles.

The scan enumerates the2^18 SOURCE x values, not2^73 relaxed parameters.
It constructs every tangent r,s for each rational fiber, checks necessary
product spans, independently replays retained literal norms, and tests full
multiplicity-safe support by R^8=0. Deduplicate r,s across repeated roots.
Returned norm candidates still require source permutations and original
Semaev fiber reconstruction before any original all-roots conclusion.

20-second native search cap and30-second supervisor,4GiB,oneCPU. A budget
pause returns visited-source counts, verified candidates as lower bounds,
null exact counts and false branch-completion flags. Complete output is
serialized before acceptance; late output is demoted and the supervisor
charges complete process lifetime. No inconclusive/timeout speed ratio.

This is a coverage derivation and native classical tangency/norm integration,
not established novelty. No external or matched F4/F5/F6 timing comparison.

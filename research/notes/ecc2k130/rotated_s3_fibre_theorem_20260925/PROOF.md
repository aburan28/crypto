# Exact finite fibres for recursive S3 on the Koblitz source curve

Let `K` be any finite field of characteristic two and
`E/K: y²+xy=x³+1`. For `a in K` put
`L(a)={P in E(K): P is affine and x(P)=a}`. The inverse of `(x,y)` is
`(x,x+y)`. Since x=0 forces y=1, `L(0)={(0,1)}`; for any nonzero a,
`L(a)` is either empty or a pair `{P,-P}`. The only affine 2-torsion
point is `(0,1)`. Use the exact polynomial

`S3(a,b,c)=(ab+ac+bc)²+abc+1 = (a+b)²c²+ab c+(ab)²+1`.

The standard summation-polynomial identity says that if affine
`P,Q,R` satisfy `P+Q+R=O`, then `S3(x(P),x(Q),x(R))=0`; equivalently,
for affine `P+Q` its x-coordinate is an S3 root because
`x(-(P+Q))=x(P+Q)`. This is the direct direction of
[Kosters–Yeo, Proposition 2.1](https://arxiv.org/pdf/1503.08001).
The reverse direction **under rational fibre membership** follows from
the complete degree count below; it does not assume rationality of a
general algebraic root.

## Local fibre theorem

Assume `L(a)` and `L(b)` are nonempty. Then the set of **distinct finite**
roots of `S3(a,b,X)` in `K` is exactly

`{x(P+Q): P in L(a), Q in L(b), P+Q != O}`.

* If `a=b=0`, both factors must be `(0,1)` and their sum is O. The
  polynomial is the constant 1, so it has no finite root.
* If `a=b!=0`, choose `P in L(a)={P,-P}`. Equal signs sum to `±2P`,
  whose common x is finite because P is not 2-torsion. Opposite signs
  sum to O. The coefficient `(a+b)²` vanishes but `ab=a²` does not,
  so S3 is linear with exactly that doubling root.
* If `a!=b` and `ab=0`, one factor is the self-inverse 2-torsion
  point. Every signed choice has the same finite sum x. S3 has
  nonzero quadratic coefficient `(a+b)²` and zero linear coefficient,
  hence is a square with precisely that one distinct root in K.
* If `a!=b` and `ab!=0`, choose `P in L(a)` and `Q in L(b)`.
  `P+Q` and `P-Q` are both affine because a and b differ. Their x's
  cannot be equal: equal x would imply `P+Q=±(P-Q)`, hence `2Q=O`
  or `2P=O`, contradicting `a,b!=0`. Both x's are S3 roots by the
  summation identity. S3 has degree two and nonzero derivative `ab`,
  so these are all its roots, both rational and rationally liftable.

This proves the four cases over every such K. In particular a
trace-one normalized quadratic cannot arise from two rationally
liftable fibres, although it can arise from arbitrary field a,b.
The implication fails without the fibre assumption: a nonliftable
`a!=0` still gives a formal linear S3 root when b=a, with no rational
factor point. It also fails for a **restricted sign set**: an S3 root
may require one of the excluded signs.

## Affine-chain theorem

Let `D_i` be factor domains formed by the **entire** rational fibres
`L(x)` for their admitted x values. Let `a_0,...,a_(m-1)` be admitted
factor x values. Suppose finite field values `c_2,...,c_(m-1)` obey

`S3(a_0,a_1,c_2)=0`,
`S3(c_j,a_j,c_(j+1))=0` for `2<=j<=m-2`, and
`S3(c_(m-1),a_(m-1),x(R))=0`

for a fixed affine rational target `R`. Then there are
`P_i in D_i` whose prefix sums of lengths `2,...,m-1` are all
affine with x-coordinates `c_2,...,c_(m-1)` and whose full sum is
**exactly R**.

To prove this, the local theorem gives signs for the first two
factors that realize `c_2`. Inductively suppose a signed tuple realizes
the prefix point `U` at `c_j`. The local theorem gives `U' in L(c_j)`
and `P_j in L(a_j)` with finite sum at `c_(j+1)`. Since
`L(c_j)={U,-U}` (possibly a singleton when `c_j=0`), either `U'=U`
or replacing **every earlier factor** by its negative changes U to
U'. All earlier prefix x values remain unchanged, and sign-complete
domains retain those factors. Append `P_j`; induction continues.
The terminal S3 step works the same way. Its finite full sum has
x-coordinate `x(R)`, hence is R or -R. If it is -R, negate **all**
factors. All prefix x values and factor memberships remain fixed,
and the full sum becomes R. When `x(R)=0`, R=-R so no choice is needed.

Conversely, any signed tuple with these affine prefix sums and affine
full sum supplies such a path by the summation identity. Thus the
affine recursive-S3 system and the exact signed-point PDP have the
same solutions after projecting to factor x tuples **provided every
relevant prefix and target is affine**. A point at infinity has no
x-coordinate. If a signed witness has an O prefix (or the target is O),
the affine system may omit that witness; an explicit identity/inverse
branch is required to claim full PDP equivalence. Other sign choices
for the same factor x tuple can still have an affine path, so an O
witness alone does not prove that its x mask is absent.

The theorem is about this S3 chain with finite K-valued intermediate
variables and complete rational sign fibres. It does not cover direct
higher S_m resultants, a denominator-cleared or incomplete Boolean
exporter, factor sets that exclude a sign, or a solver's treatment of
timeouts. In particular it supplies no solver complexity or ECDLP
speed claim. The exhaustive small-field controls in this PR test the
finite cases and the induction against a separately implemented group
law; they are a check on the derivation, not a substitute for it.

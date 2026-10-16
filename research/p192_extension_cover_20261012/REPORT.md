# Exact extension-cover screen for standard NIST P-192

**Latest relation criterion:** The [signed-base scalar-anchor audit](reports/run-20261015T1857Z/REPORT.md)
enumerates all 695,520 three-sums in the fixed 160-point Type I quotient
base. It finds no overlap between independently known-scalar anchors and
nonzero-unknown-row targets, and every repeated rank-bearing target has the
same full formal signed-point vector. The criterion and exact receipt are
linked there.

**Adaptive-base continuation:** The [64-query rank-balance audit](reports/run-20261015T1926Z/REPORT.md)
constructs one new quotient point per known-scalar equation and checks every
point against the old three-sum support. The 64 forced rows bring 64 private
columns and leave the 16 old logarithms unconstrained. Its Lean theorem
formalizes the unique extension of any old-log assignment.

**Short-sum collision screen:** The [new exact audit](reports/run-20261015T1955Z/REPORT.md)
tests the 64 selected adaptive points against all two- and four-sums of the
160 signed old classes. It finds zero hits in 64 direct pair lookups and
824,320 pair-complement tests, independently replayed in Python. Together
with the earlier one- and three-sum tests, the minimum old-base sum length
for each selected new point is at least five. Lean checks the group identity
that would turn such a collision into an old-base equation.

**Divisor-support geometry:** The [split-support audit](reports/run-20261015T2042Z/REPORT.md)
checks all 160 fixed-origin degree-three residual divisors. Their
\(x\)-polynomials are pairwise coprime across 12,720 class pairs; 20 split
fully into rational points. Restricting the collector to those 20 classes
leaves 143 geometry-accepted rank-bearing three-sum targets, with zero hits
in the frozen 64 fresh known-scalar queries. Python and PARI independently
replay the support and quotient-group counts, and Lean checks the support
disjointness implication and exact uniform-query inequality.

**Continuation:** The [three-branch follow-up](BRANCH_CONTINUATION.md) proves
that the nonsplit-quadratic direct genus-3 branch is excluded as well,
constructs a descended Type I nonhyperelliptic genus-3 cover and a
nonfixed-coefficient Scholten genus-2 descent, and verifies same-field
genus-3 and genus-2 controls. The original screen and
receipts below remain the baseline for those results.
The continuation also gives a prime-field Type I discriminant criterion and
constructs exact P-224 and P-256 subgroup conorms whose three-step Frobenius
norms are nonzero over \(\mathbf F_{p^2}\).

**Result.** For standard P-192 viewed over \(\mathbf F_{p^6}\), every Legendre
cross-ratio of its nonzero 2-torsion has norm \(-1\) to \(\mathbf F_{p^2}\).
Thus it cannot be written in the split-linear-root hyperelliptic genus-3
family \(y^2=x(x-\alpha)(x-\alpha^{p^2})\) used in the explicit cubic-cover
construction of Joux and Vitse. A separate, exact genus-2 cover
\(H:v^2=(u^2)^3+a u^2+b\), \((u,v)\mapsto(u^2,v)\), does preserve the P-192
prime-order subgroup, but is already defined over \(\mathbf F_p\). The
subsequent [branch continuation](BRANCH_CONTINUATION.md) excludes the direct
quadratic-root genus-3 family, constructs the nonhyperelliptic Type I descent,
and verifies two further subgroup-preserving controls. It also constructs a
nonfixed-coefficient Scholten descent. Both descended Jacobian relation
pipelines remain open.

## Instance and exact receipts

The instance is NIST P-192 / secp192r1, with

```
p = FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
a = p - 3
b = 64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1
r = FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831
Gx = 188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012
Gy = 07192B95FFC8DA78631011ED6B24CDD573F977A11E794811
```

The checked Sage run independently computed \(\#E(\mathbf F_p)=r\), proved
\(r\) prime, verified \(rG=0\), proved the 2-division cubic irreducible,
and checked all six Legendre norms. It constructed the genus-2 Jacobian
divisor \(D=[u^2-G_x,G_y]\), checked \(D\ne0\) and \(rD=0\), and used the
Frobenius trace recurrence to compute exact extension orders. The measured
14.297 s is a diagnostic duration for this screen, **not** an attack time.
The accepted Sage runtime and module receipt is saved in
[`sage-runtime-info.txt`](sage-runtime-info.txt); the machine-readable result
is [`results.json`](results.json).

| Extension degree \(k\) | Field bits | \(v_2(\#E(\mathbf F_{p^k}))\) | \(r\mid p^k-1\) |
| ---: | ---: | ---: | :---: |
| 1 | 192 | 0 | no |
| 2 | 384 | 0 | no |
| 3 | 576 | 2 | no |
| 5 | 960 | 0 | no |
| 6 | 1152 | 4 | no |

The last column screens the direct MOV/Frey–Rück pairing field at these
specific degrees; it says nothing about larger embedding degrees. The
point-count recurrence is \(T_0=2\), \(T_1=p+1-r\),
\(T_k=T_1T_{k-1}-pT_{k-2}\), and \(N_k=p^k+1-T_k\).

## Norm obstruction theorem

Let \(p\) be odd and let \(E/\mathbf F_p\) have irreducible 2-division cubic
with distinct roots \(e_0,e_1,e_2\) in \(\mathbf F_{p^3}\), numbered so
\(e_i^p=e_{i+1}\) cyclically. Put \(q=p^2\) and
\(\lambda=(e_2-e_0)/(e_1-e_0)\). Then

\[
N_{\mathbf F_{p^6}/\mathbf F_{p^2}}(\lambda)
=\lambda\lambda^q\lambda^{q^2}=-1.
\]

The other five cross-ratios obtained by permuting the roots also have norm
\(-1\). Consequently this elliptic curve, after base change to
\(\mathbf F_{p^6}\), has no Legendre model in the split-linear-root family
\(y^2=x(x-\alpha)(x-\alpha^q)\): its corresponding ratio
\(\alpha^q/\alpha\) would have norm \(1\).

**Proof.** The \(q\)-Frobenius reverses the cyclic ordering of the three
roots. The three conjugates are

\[
\frac{e_2-e_0}{e_1-e_0},\qquad
\frac{e_1-e_2}{e_0-e_2},\qquad
\frac{e_0-e_1}{e_2-e_1}.
\]

Their product is \(-1\) by cancellation. Any permutation gives the same
cyclic product with the labels renamed. A split-linear-root model has roots
\(0,\alpha,\alpha^q\), hence a cross-ratio \(\alpha^q/\alpha\), whose norm
is \(\alpha^{q^3}/\alpha=1\). The cross-ratio is unchanged by the affine
change of the 2-division roots that accompanies a Weierstrass isomorphism.
For P-192, the exact irreducibility check supplies the hypotheses.

This excludes the **split-linear-root subfamily**. In the terminology of
[Joux–Vitse, §4.1 and Appendix A](https://eprint.iacr.org/2011/020.pdf),
their general hyperelliptic genus-3 equation also permits a quadratic
\(h(x)\); the theorem does not classify that branch. Their §4.1 also treats
nonhyperelliptic genus-3 covers, and §4.2 treats Scholten genus-2 covers.
Those remain separate gates.

## Exact subgroup-preserving genus-2 control

Let \(F(x)=x^3+ax+b\). Because \(b\ne0\), \(F\) is squarefree, and \(p\)
is odd, \(F(u^2)\) is squarefree of degree six. Thus
\(H:v^2=F(u^2)\) is smooth of genus two. The map
\(\pi:H\to E\), \((u,v)\mapsto(u^2,v)\), has degree two over
\(\mathbf F_p\). For the P-192 generator, the Mumford class
\(D=[u^2-G_x,G_y]=\pi^*(G-O)\) is nonzero and has order \(r\), as checked
by Sage. The divisor identity \(\pi_*\pi^*=[2]\) gives
\(\pi_*D=2G\). Since \(r\) is odd, pullback is injective on the
\(r\)-subgroup. After base change to \(\mathbf F_{p^6}\), the
\(\mathbf F_{p^6}/\mathbf F_{p^3}\) norm of this
\(\mathbf F_p\)-rational class is \(2D\), and
\(\pi_*N(D)=4G\).

This validates the geometry and subgroup-map gate. It provides a control
for candidate-cover software: a mere lower-field model or nonzero subgroup
image must be distinguished from a new transfer route. Here the source
curve, cover and divisor all descend to the **original** prime field.

## Formal proof map

The pinned Lean 4 project is [`lean/`](lean/). Reproduce with
`cd research/p192_extension_cover_20261012/lean && lake build CrossRatio`.
The successful build and `#print axioms` output are in
[`lean/lean-build.log`](lean/lean-build.log). No `sorry`, `admit`, custom
axiom, or native evaluation appears in the theorem module.

| Mathematical statement | Lean declaration | Status and exact boundary |
| --- | --- | --- |
| Product of cyclic cross-ratios is \(-1\) in any field | `P192ExtensionCover.cyclic_cross_ratio_product` | Checked; `propext`, `Classical.choice`, `Quot.sound` only. |
| Product cannot equal one in characteristic other than two | `P192ExtensionCover.cyclic_cross_ratio_not_norm_one` | Checked with explicit hypothesis \(-1\ne1\); same foundations. |
| Type I discriminant factor is nonzero when minus one is nonsquare | `P192ExtensionCover.type_i_discriminant_factor_ne_zero` | Checked under explicit field hypotheses; finite-field and cover criteria are detailed in the continuation. |
| Push-pull equal to doubling implies injective pullback where doubling is injective | `P192ExtensionCover.pullback_injective_of_push_pull_two` | Checked; no axioms. |
| Full P-192 finite-field and cover theorem | Open | Formalize the P-192 field, Frobenius permutation, norm identification, smooth curve/map, and divisor push-pull identity. Sage verifies the exact instance; these bridges are not yet Lean declarations. |

The Lean theorem is the algebraic core of the norm obstruction, not a
formal certificate for the complete elliptic-curve statement. The written
proof above gives the structural steps; the last row identifies the precise
formalization work that remains.

## Cost and next construction gates

The explicit Joux–Vitse factor base in their cubic-cover experiment has
size on the order of \(p\). Applied unchanged at P-192 scale, even
enumerating that base would require on the order of \(2^{192}\) entries,
while a generic single-target reference has square-root scale near
\(2^{96}\) group operations. These are different operation units, so the
numbers are a feasibility screen for that **particular base policy**, not a
measured speed ratio or a lower bound on every cover method. See
[Joux–Vitse, §§3–4](https://eprint.iacr.org/2011/020.pdf) for the
decomposition method and base choice.

The [branch continuation](BRANCH_CONTINUATION.md) resolves the direct
quadratic-root test and constructs Type I genus-3 and nonfixed Scholten
genus-2 covers over their lower fields. It gives explicit lower-field Mumford
coordinates for a nonzero P-192 Scholten subgroup class, 16 verified
subgroup-usable classes in a compact base, and a 64-query zero-yield relation
probe with an exact three-sum coverage bound. A second known-scalar policy
adds tracked base-class shifts and yields 0 relations in 64 queries, rank
0/16; 64 fresh geometric cover classes also yield 0 exact three-sum matches.
The latter class logs are unknown, so those matches would be decomposition
observations rather than relation equations. Lean checks the shifted relation
identity and its unchanged 816-image support bound. The continuation also gives
checked Type I generator fibers and a conorm class of exact P-192 subgroup order, proved
through degree-two push-pull and checked quotient images. Three elliptic
quotients reduce its genus-three subgroup class to
\((-2G,0,0)\); the three-step Frobenius norm reduces to
\((-2G,-2G,-2G)\), giving a nonzero class of exact order \(r\) over
\(\mathbf F_{p^2}\). A new frozen Type I probe builds 16 distinct nonzero
diagonal \(r\)-torsion vectors from \(\mathbf F_p\) cover points in these
three quotient coordinates. It finds 0 exact three-sum matches and quotient
incidence rank 0/16 in 64 known-scalar queries; an independent verifier
enumerates all 816 triples. Lean checks the six-periodic Frobenius eigenvector
step that explains diagonality when \(p^6\ne1\pmod r\); the finite-field and
Jacobian bridges are stated separately. P-224 and P-256 have the same checked diagonal norm
pattern, using \(G\) and \(17G\), respectively; independent point-coordinate
verifiers replay the quotient values. A separate
\(\mathbf F_{23^3}\) Scholten instance verifies the complete relation,
rank, target, and same-point rho pipeline. P-192 final rank, target descent,
and matched one-target cost remain unknown for the tested compact policy;
the Type I branch still needs a relation method and complete target cost;
generic plane-quartic divisor reduction remains available as a separate
implementation route.

The [fixed-origin nonspeciality certificate](reports/run-20261015T2117Z/REPORT.md)
now checks all 160 signed Type I degree-three representatives. Their
squarefree residual (x)-cubics and nonzero quadratic (y)-coefficients
make every residual triple noncollinear, so Riemann--Roch gives a unique
effective degree-three divisor in each fixed-origin class. This identifies
the exact geometry constraint on a collector that tries to vary those
representatives while keeping the same origin; it does not replace the
separate generic Jacobian arithmetic gate.

The [line-section continuation](reports/run-20261015T2201Z/REPORT.md)
enumerates all 9,870 pair-lines through the 141 rational support points,
finds 9,784 fresh rational intersections, and independently projects all
9,925 selected points to distinct first-quotient subgroup classes. A
non-splitting line through two fixed support points, together with the exact
quotient-base cover anchor, gives nonspecial fixed-origin degree-three source
representatives for all selected points. On the 141-class projected base,
the frozen uniform and character-conditioned 64-query known-scalar laws
each yield zero three-sum matches and rank zero. The first 16 line-derived
classes have 816 distinct three-sum targets, 206 admitted by the two-fiber
character law. The [all-character quotient reduction and collector run](reports/run-20261015T2235Z/REPORT.md)
reduces all 141 cofactor-projected quotient classes to explicit nonspecial
fixed-origin degree-three Jacobian representatives. The exact first-16
three-sum support splits as 206, 218, 199, and 193 across the four fiber
character cells. A balanced four-cell known-scalar policy admits all 816
support targets and fills 64 queries after 76 scalar trials, compared with
241 for a ++-only control. Both 64-query arms have zero full-base matches
and rank zero under the frozen 141-class pair-index lookup.

The [complete 141-class three-sum census](reports/run-20261015T2255Z/REPORT.md)
finds 10,011 distinct pair targets and 477,191 distinct nonidentity
three-sum targets, saturating the unordered-three-sum capacity.
Independent integer-field replay matches the full support digest.
For an ideal uniform nonzero scalar query, the exact overlap probability
is \(477191/(r-1)\); conditional character-cell hit rates require
cell-population denominators in addition to eligible support counts.
The general capacity bound is formalized in Lean. A frozen 512+512
support-versus-uniform sample estimates the ++ cell hit-rate multiplier
at 1.048 and the balanced four-cell multiplier at about 1.002; these
exploratory estimates do not resolve a useful conditional enrichment.

The [signed short-support screen](reports/run-20261015T2316Z/REPORT.md)
extends that census to all positive and negated supports of lengths one,
two, and three. Its 12 zero-target and support-intersection cells are
empty, so no signed relation row of this declared shape arises on the
141-class base. Both integer-field and Sage group laws replay the
support sets; longer supports and known-scalar target overlap remain
separate collector measurements.

The [split-line leaf-rank analysis](reports/run-20261015T2346Z/REPORT.md)
uses the 4,892 split line cuts to expose 9,784 additional distinct
quotient points. Each cut has two private fresh-point columns, so the
geometric line equations are independent but no nonzero combination
eliminates all fresh columns to constrain the old 141 alone. The first
32 constant-sum equations pass independent Sage and prime-field checks;
the full incidence and distinctness are frozen. Even at the combined
9,925-point size, the unordered three-sum ceiling is 162,993,977,475,
giving a uniform 64-query union bound of approximately (2^{-148.754}).
The private-column implication has a checked Lean proof.

The [second-generation line-closure census](reports/run-20261015T2357Z/REPORT.md)
tests all 11,040 pair lines containing one of the first 64 new points
against the full 9,925-point coordinate set. It finds 160 closed cuts,
all alternate presentations of the 32 original quartets, and zero new
closed quartets. Independent GP and prime-field implementations agree;
a fresh-checked Lean theorem shows duplicate rows cannot enlarge the
relation span.

The [full split-line constant replay](reports/run-20261016T0116Z/REPORT.md)
projects all 9,784 new points and verifies that every one of the 4,892
quartets sums to the same first-quotient subgroup point \(C\). Its unique
fourth part \(D\) in the odd-order subgroup centers each row to a zero
sum without changing point-label distinctness or the private-column
obstruction. Independent prime-field arithmetic replays every sum and
the full coordinate digest. A checked Lean theorem proves the abstract
four-point centering identity; the exact finite-field quotient model
remains a concrete Sage and independent-integer replay. Since centering
a known-scalar target by \(D\) also changes its scalar anchor, a
collector must account explicitly for the scalar of \(D\) or carry its
unknown as another column. A cached quartic inverse matches all three
quotient images on 80 frozen points; its local timing is diagnostic
pending isolated paired measurement.

The [complete 9,925-point collinearity census](reports/run-20261016T0141Z/REPORT.md)
visits all 49,247,850 unordered pairs of frozen affine quartic points.
Exactly the 4,892 original split quartets contain three or more selected
points; every one has four, with two old and two line-derived points.
An independent reverse-order, block-inversion implementation reproduces
every canonical line key and incidence set. The full selected-endpoint
tangency census finds zero tangent pairs. For any other line with four
distinct intersection points, at least two lie outside the frozen set;
a fresh-checked Lean theorem proves this implication from the complete
census and distinct-cut hypotheses. Residual ramification and reuse of
points introduced by later cuts remain separate collector questions.

## Reproduction

From this repository's root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/sage-runtime-info.txt
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/screen.py > research/p192_extension_cover_20261012/sage-screen.log 2>&1
cd research/p192_extension_cover_20261012/lean
lake build CrossRatio > lean-build.log 2>&1
```

The exact parameters above are from [NIST FIPS 186-4, Appendix D.1.2.1](https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf).

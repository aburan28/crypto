# Exact extension-cover screen for standard NIST P-192

**Continuation:** The [three-branch follow-up](BRANCH_CONTINUATION.md) proves
that the nonsplit-quadratic direct genus-3 branch is excluded as well,
constructs a descended Type I nonhyperelliptic genus-3 cover and a
nonfixed-coefficient Scholten genus-2 descent, and verifies same-field
genus-3 and genus-2 controls. The original screen and
receipts below remain the baseline for those results.

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
quadratic-root test, constructs a Type I genus-3 cover and a nonfixed
Scholten genus-2 cover over their lower fields, and verifies the earlier
same-field controls. The next gate is an explicit Jacobian divisor image and
conorm-norm replay on each descended cover. Their subgroup-usable factor bases,
relation yields, ranks, final linear algebra, target descent, and matched single-target rho reference
remain **unknown** for the P-192 instance.

## Reproduction

From this repository's root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/sage-runtime-info.txt
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/screen.py > research/p192_extension_cover_20261012/sage-screen.log 2>&1
cd research/p192_extension_cover_20261012/lean
lake build CrossRatio > lean-build.log 2>&1
```

The exact parameters above are from [NIST FIPS 186-4, Appendix D.1.2.1](https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf).

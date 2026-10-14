# Three P-192 cover branches: exact direct-family gate and verified controls

**Result.** Standard NIST P-192 has no direct Joux–Vitse hyperelliptic genus-3 model over the tower \(\mathbb F_{p^6}/\mathbb F_{p^2}/\mathbb F_p\), including the nonsplit-quadratic-root branch. The earlier cross-ratio norm excludes the linear and split-quadratic cases; the nonsplit branch has \(j\notin\mathbb F_p\), while P-192 has \(j\in\mathbb F_p\). We also constructed a smooth, nonhyperelliptic genus-3 plane quartic over \(\mathbb F_p\) with a degree-two map to P-192, and an explicit Scholten-form genus-2 cover over \(\mathbb F_{p^3}\). The latter transfers the P-192 generator to a verified nonzero order-\(r\) Jacobian class. These two constructions are controls over the original curve's field or a larger one; neither is yet the genuinely descended Joux–Vitse cover-and-decomposition route.

The exact inputs are [in the earlier report](REPORT.md), and the executable receipt is [`branch_gates.py`](branch_gates.py) with [JSON](branch_gates.json), [Sage runtime](branch_sage_runtime_info.txt), and [run log](branch_gates.log). The subgroup is the standard 192-bit prime-order P-192 group; \(p=6277101735386680763835789423207666416083908700390324961279\), \(r=6277101735386680763835789423176059013767194773182842284081\).

## Complete direct hyperelliptic genus-3 exclusion

Let \(p>3\), \(q=p^2\), \(K=\mathbb F_{p^6}\), and let \(d\) be a nonsquare in \(\mathbb F_q\). For a nonsingular quartic in the nonsplit branch
\[
 C_\alpha:y^2=(x^2-d)(x-\alpha)(x-\alpha^q),\qquad \alpha\in K\setminus\mathbb F_q,
\]
its elliptic \(j\)-invariant does not lie in \(\mathbb F_p\). This strengthens the previous P-192 screen from one split-linear-root subfamily to the complete **direct** hyperelliptic genus-3 family of [Joux–Vitse, §4.1 and Appendix A](https://iacr.org/archive/eurocrypt2012/72370010/72370010.pdf). It does not decide whether a different curve in P-192's \(\mathbb F_{p^6}\) isogeny class has such a model.

**Proof.** Choose \(\beta^2=d\) in \(L=\mathbb F_{p^{12}}\). The branch-root cross-ratio gives \(\lambda\in\mu_N\setminus\{1\}\), where \(N=p^4-p^2+1\), and
\(j=256(\lambda^2-\lambda+1)^3/[\lambda^2(\lambda-1)^2]\). These identities and the constructive Cayley parameterization are proved in the existing [two-branch theorem](../iso1_weak_classes_20261007/two_branch_20261009/REPORT.tex). Since \(N\mid p^6+1\), \(\lambda^{p^6}=\lambda^{-1}\). Among the six Legendre transforms with equal \(j\), only \(\lambda\) and \(\lambda^{-1}\) remain in \(\mu_N\): every other transform has norm to \(K\) equal to \(2-(\lambda+\lambda^{-1})\) or its inverse. Norm one would force \(\lambda^2-\lambda+1=0\), hence \(\lambda^6=1\). This is impossible because \(N\) is odd and \(N\equiv1\pmod3\). If \(j\in\mathbb F_p\), then \(j(\lambda^p)=j(\lambda)\), so \(\lambda^p=\lambda\) or \(\lambda^{-1}\). Yet \(N\equiv1\pmod{p-1}\) and \(N\equiv1\pmod{p+1}\); either equality forces \(\lambda=1\), a contradiction. The cleared six-factor identity is kernel-checked in Lean as `P192ExtensionCover.legendre_j_difference_factor`. The remaining finite-field torus and norm bridge is written here and in the prior theorem supplement, but is not yet a full Lean theorem.

A linear \(h\), or split squarefree quadratic \(h\), reduces by an \(\mathbb F_q\) coordinate change to the cubic family in [the earlier report](REPORT.md); its required norm-one cross-ratio conflicts with P-192's norm \(-1\). Thus the direct family is exhausted. The exact P-192 \(j\)-invariant is \(6234286251230310114240839169629130138801351179850969208331\) in \(\mathbb F_p\); the run also checked it differs from 0 and 1728.

## Explicit nonhyperelliptic genus-3 cover

For \(E:y^2=x^3+ax+b\), define the projective plane quartic
\[
 C: Z^4=X^3W+aXW^3+bW^4,
 \qquad \pi(x,z)=(x,z^2).
\]
It is smooth because, in the affine chart, simultaneous vanishing of its derivatives forces \(z=0\) and a common root of \(f(x)=x^3+ax+b\) and \(f'(x)\); the previous Sage receipt checked \(f\) squarefree. At infinity its sole point \([1:0:0]\) has nonzero \(W\)-derivative. A smooth plane quartic has genus 3 and is nonhyperelliptic. The function-field extension adjoining a square root of \(y\) has degree two. For a degree-two cover, \(\pi_*\pi^*=[2]\), so its pullback is injective on P-192's odd-order \(r\)-subgroup. The polynomial map identity is checked in Lean as `P192ExtensionCover.quartic_projection_identity`.

This construction is defined already over \(\mathbb F_p\). Its existence establishes a nonhyperelliptic cover and subgroup transfer, but does not establish the lower-field \(\mathbb F_{p^2}\) cover classified in the Joux–Vitse/Momose–Chao branch. A literal plane-quartic Jacobian class for \(G\), a \(p^2\)-field descended model, and a Diem relation solver remain to be built.

## Explicit Scholten-form genus-2 cover

Let \(\zeta\) be a root of \(f\) in \(\mathbb F_{p^3}\), and put
\[
 s=-\frac{3\zeta^2+a}{3\zeta},\quad A=s^3,\quad B=3s^3+3\zeta s^2.
\]
Here \(\zeta\ne0\) and \(f'(\zeta)\ne0\), so \(s\ne0\). The checked identity
\[
 f(\zeta+s(t+1))=At^3+Bt^2+Bt+A
\]
puts P-192 in the Scholten form \(y^2=At^3+Bt^2+\sigma(B)t+\sigma(A)\) for \(\sigma:x\mapsto x^{p^3}\), since \(A,B\in\mathbb F_{p^3}\). Define
\[
 H:v^2=Au^6+Bu^4+Bu^2+A,
 \qquad (u,v)\longmapsto(\zeta+s(u^2+1),v)\in E.
\]
Sage checked the sextic squarefree and genus two. It constructed \(D=[u^2-t_G,G_y]\) in \(\operatorname{Jac}(H)(\mathbb F_{p^3})\), where \(t_G=(G_x-\zeta)/s-1\), and verified \(D\ne0\) and \([r]D=0\). The map has degree two and \(\pi_*D=2G\). Over \(\mathbb F_{p^6}/\mathbb F_{p^3}\), the norm of this already rational class is \(2D\); pushforward gives \(4G\). The palindromic identity is kernel-checked in Lean as `P192ExtensionCover.scholten_palindromic_identity`.

This is the fixed-coefficient Scholten specialization: the cover and map are defined over \(\mathbb F_{p^3}\), and P-192 itself is defined over \(\mathbb F_p\). The explicit coefficients in the polynomial basis \(1,\zeta,\zeta^2\) are in the JSON receipt. The route demonstrates literal geometry and subgroup preservation; the nontrivial coefficient-twisted Scholten descent and its effect on factor-base relations remain open.

## Measured gates and complete-cost boundary

The frozen 256-abscissa sample uses Python MT19937 `randrange(p)` with seed `0x19220261013`; its exact decimal-list SHA-256 is in the JSON receipt. The Sage run completed in 4.487 s on this shared host. These timings are diagnostics, with no isolated CPU timing ratio.

| Branch | Exact direct geometry | Transferred \(r\)-class | 256 abscissae with candidate points | Relation yield / rank / LA / target descent | Full one-target cost |
| --- | --- | --- | ---: | --- | --- |
| Nonsplit quadratic genus 3 | Excluded for P-192 by \(j\)-degree | Not applicable | Not applicable | Not applicable | Not applicable for the direct model |
| Nonhyperelliptic genus 3 | Smooth degree-two plane quartic over \(\mathbb F_p\) | Follows from push-pull; explicit Jacobian arithmetic open | 115; 0.003 s eligibility check | Unmeasured | Unknown |
| Scholten genus 2 | Smooth degree-two cover over \(\mathbb F_{p^3}\) | \(D\ne0, [r]D=0\) checked | 139; 0.005 s eligibility check | Unmeasured | Unknown |

The sample tests only square/fourth-power eligibility of chosen \(x\)-coordinates. It does not enumerate a distinct usable factor base, collect a relation, solve a matrix, or descend an unseen target. The original Joux–Vitse factor-base policy has \(x\in\mathbb F_p\), so a direct complete scan at P-192 has \(p\) candidate abscissae, about \(2^{192}\). Its reported practical run instead used a roughly 26-bit base prime and a 149-bit subgroup over a degree-six extension [Joux–Vitse, §§3–5](https://iacr.org/archive/eurocrypt2012/72370010/72370010.pdf). The P-192 full solve is therefore not a completed measurement; for either surviving control the target-dependent online time and paired same-point rho time remain unknown under the repository's one-target measurement contract.

## Formal proof map and reproduction

`lean/CrossRatio.lean` builds under Lean 4.34.1 and pinned mathlib revision `d13f23b723b8a846827a245b89c10fc7d3f11612`. Its `#print axioms` output is in `lean/lean-build.log`; each new theorem uses only `propext` and/or `Quot.sound`. The new polynomial identity, Scholten coordinate identity, and quartic projection identity are fully checked. The finite-field norm and cross-ratio classification, P-192 elliptic isomorphism, smooth plane-curve genus, divisor push-pull, and Sage Jacobian arithmetic are outside those Lean declarations. A full Lean theorem for the direct-family exclusion and either cover is still open.

Run from the `crypto` worktree root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/branch_sage_runtime_info.txt
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/branch_gates.py > research/p192_extension_cover_20261012/branch_gates.log 2>&1
(cd research/p192_extension_cover_20261012/lean && lake build CrossRatio > lean-build.log 2>&1)
```

The next useful implementation is a literal \(\mathbb F_{p^2}\) nonhyperelliptic descended cover and a nontrivial \(\mathbb F_{p^3}\) Scholten twist, each with explicit conorm-norm replay. Those constructions must precede any relation solver or full cost. A bounded small-field end-to-end control can then exercise factor-base construction, ordinary relation yield, rank, final linear algebra, and target descent before any P-192 extrapolation.

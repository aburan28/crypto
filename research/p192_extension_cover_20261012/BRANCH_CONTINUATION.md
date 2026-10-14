# Three P-192 cover branches: exact direct-family gate and verified controls

**Result.** Standard NIST P-192 has no direct Joux–Vitse hyperelliptic genus-3 model over the tower \(\mathbb F_{p^6}/\mathbb F_{p^2}/\mathbb F_p\), including the nonsplit-quadratic-root branch. The earlier cross-ratio norm excludes the linear and split-quadratic cases; the nonsplit branch has \(j\notin\mathbb F_p\), while P-192 has \(j\in\mathbb F_p\). We also constructed a smooth, nonhyperelliptic genus-3 plane quartic over \(\mathbb F_p\) with a degree-two map to P-192, and an explicit Scholten-form genus-2 cover over \(\mathbb F_{p^3}\). The latter transfers the P-192 generator to a verified nonzero order-\(r\) Jacobian class. These two constructions are controls over the original curve's field or a larger one; neither is yet the genuinely descended Joux–Vitse cover-and-decomposition route.

The exact inputs are [in the earlier report](REPORT.md), and the executable receipt is [`branch_gates.py`](branch_gates.py) with [JSON](branch_gates.json), [Sage runtime](branch_sage_runtime_info.txt), and [run log](branch_gates.log). The subgroup is the standard 192-bit prime-order P-192 group; \(p=6277101735386680763835789423207666416083908700390324961279\), \(r=6277101735386680763835789423176059013767194773182842284081\).

The distinction between a subgroup-preserving map and an effective descent is central here. [Momose–Chao, §1](https://eprint.iacr.org/2009/236.pdf) formulates the favorable cover case with an induced isogeny from the cover Jacobian to the Weil restriction of the elliptic curve. The two explicit controls below verify preservation of the selected rational P-192 subgroup, but no such Weil-restriction isogeny or lower-field relation algorithm has been verified for them. Their existence is therefore a geometry and arithmetic control for the requested cover search.

## Complete direct hyperelliptic genus-3 exclusion

Let \(p>3\), \(q=p^2\), \(K=\mathbb F_{p^6}\), and let \(d\) be a nonsquare in \(\mathbb F_q\). For a nonsingular quartic in the nonsplit branch
\[
 C_\alpha:y^2=(x^2-d)(x-\alpha)(x-\alpha^q),\qquad \alpha\in K\setminus\mathbb F_q,
\]
its elliptic \(j\)-invariant does not lie in \(\mathbb F_p\). This strengthens the previous P-192 screen from one split-linear-root subfamily to the complete **direct** hyperelliptic genus-3 family of [Joux–Vitse, §4.1 and Appendix A](https://iacr.org/archive/eurocrypt2012/72370010/72370010.pdf). It does not decide whether a different curve in P-192's \(\mathbb F_{p^6}\) isogeny class has such a model.

**Proof.** Choose \(\beta^2=d\) in \(L=\mathbb F_{p^{12}}\). The branch-root cross-ratio gives \(\lambda\in\mu_N\setminus\{1\}\), where \(N=p^4-p^2+1\), and
\(j=256(\lambda^2-\lambda+1)^3/[\lambda^2(\lambda-1)^2]\). These identities and the constructive Cayley parameterization are proved in the existing [two-branch theorem](../iso1_weak_classes_20261007/two_branch_20261009/REPORT.tex). Since \(N\mid p^6+1\), \(\lambda^{p^6}=\lambda^{-1}\). Among the six Legendre transforms with equal \(j\), only \(\lambda\) and \(\lambda^{-1}\) remain in \(\mu_N\): every other transform has norm to \(K\) equal to \(2-(\lambda+\lambda^{-1})\) or its inverse. Norm one would force \(\lambda^2-\lambda+1=0\), hence \(\lambda^6=1\). This is impossible because \(N\) is odd and \(N\equiv1\pmod3\). If \(j\in\mathbb F_{p^d}\), then \(j(\lambda^{p^d})=j(\lambda)\), so \(\lambda^{p^d}=\lambda\) or \(\lambda^{-1}\). For \(d=1,2\), \(\gcd(N,p^d-1)=\gcd(N,p^d+1)=1\); for \(d=2\), the only possible common divisor of \(N\) and \(p^2+1\) is 3, and \(N\equiv1\pmod3\). Hence either equality forces \(\lambda=1\). For \(d=3\), either equality implies \(\lambda^{p^6}=\lambda\), contradicting \(\lambda^{p^6}=\lambda^{-1}\) unless \(\lambda=1\). Thus the degree of \(j\) over \(\mathbb F_p\) is exactly six. The cleared six-factor identity is kernel-checked in Lean as `P192ExtensionCover.legendre_j_difference_factor`. The remaining finite-field torus and norm bridge is written here and in the prior theorem supplement, but is not yet a full Lean theorem.

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

The exact scan count for this policy is \(p=6277101735386680763835789423207666416083908700390324961279\). Even a one-bit eligibility bitmap indexed by every candidate would require at least
\(\lceil p/8\rceil=784637716923335095479473677900958302010488587548790620160\) bytes, nearly \(2^{189}\) bytes. The candidate-inspection count exceeds \(\sqrt r\) by approximately \(2^{96}\). These are exact work and storage bounds **conditional on that full-abscissa policy**, not measured wall time or a bound on all possible factor-base constructions. The prior 256-coordinate sample estimates only eligibility; it cannot be used to infer relation yield or the number of subgroup-usable base classes.

## End-to-end Jacobian pipeline control

An executable [small-field control](small_field_pipeline.py) exercises ordinary relation collection, useful rank, final modular linear algebra, target descent, and a same-point rho reference on an explicit genus-2 Jacobian. It uses \(E/\mathbb F_{101}:y^2=x^3+x+46\) of prime order 109 and the degree-two cover \(H:v^2=u^6+u^2+46\). Sage verifies \(\#\operatorname{Jac}(H)(\mathbb F_{101})=11663=107\cdot109\), and constructs a nonzero order-109 class \(D=[u^2-2,37]\). Each factor-base class \(P_i\) comes from an affine rational point of \(H\), projected to the order-109 subgroup as \(Q_i=107P_i\). The base contains 32 distinct projected classes. For each ordinary query \(aD\), a pair index searches for an exact three-point Jacobian decomposition \(aD=P_i+P_j+P_k\); multiplying by 107 gives the modular log relation. The code verifies every divisor equality before inserting a row.

| Field and branch | Base classes | Ordinary queries | Verified relations | Independent rows / rank | Fresh-target descent | Same-point rho | Timing status |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| \(\mathbb F_{101}\), same-field genus 2 | 32 | 165 | 63 | 32 / 32 | 7 offset queries; scalar 18 replayed | scalar 18 replayed; 10 walk steps | Shared-host diagnostic |

The run records the deterministic input law, source hash, base digest, each rank-increasing relation witness, phase durations, and target point digest in [JSON](small_field_pipeline.json), with the [Sage runtime receipt](small_field_runtime_info.json), [raw log](small_field_pipeline.log), and [SHA-256 list](small_field_receipts.sha256). The target is generated only after the relation matrix and base logs are ready; its scalar is excluded from the online algorithm and checked after recovery. The online target interval was 0.837 s; the matched rho interval was 0.000391 s on this shared host. These are diagnostic timings for a 109-element group and do not constitute an isolated speed comparison. The complete P-192 pipeline remains unmeasured because its genuinely descended cover and a usable factor base have not been constructed.

An initial control run failed its scalar replay: target descent omitted multiplication by \(107^{-1}\bmod109\) after cofactor projection. The [failed receipt](small_field_pipeline_attempt1.log) is retained; the corrected formula and the group-projection and scalar identities are checked in Lean as `projected_three_point_relation` and `projected_target_scalar`. This failure illustrates why successful decomposition and full matrix rank alone are not sufficient correctness gates.

## Formal proof map and reproduction

`lean/CrossRatio.lean` builds under Lean 4.34.1 and pinned mathlib revision `d13f23b723b8a846827a245b89c10fc7d3f11612`. Its `#print axioms` output is in `lean/lean-build.log`; the declarations use standard Lean foundations only. The polynomial identity, Scholten coordinate identity, quartic projection identity, integer-scalar group projection, and field-scalar target formula are checked. The latter two formalize the correction exposed by the failed control run; Sage checks their Jacobian instantiation and scalar replay. The finite-field norm and cross-ratio classification, P-192 elliptic isomorphism, smooth plane-curve genus, divisor push-pull, and Sage Jacobian arithmetic are outside those Lean declarations. A full Lean theorem for the direct-family exclusion and either cover is still open.

| New statement | Lean declaration | Verification boundary |
| --- | --- | --- |
| Project an exact three-point relation by an integer cofactor | `P192ExtensionCover.projected_three_point_relation` | Kernel-checked for any additive commutative group; Jacobian relations separately checked by Sage. |
| Recover the scalar after cofactor projection | `P192ExtensionCover.projected_target_scalar` | Kernel-checked over any field with nonzero cofactor; the modular 109 instance and scalar replay checked by Sage. |
| Exact P-192 direct-family exclusion and effective descent | Open | The finite-field torus, elliptic-cover geometry, and Weil-restriction isogeny are not formalized in Lean. |

Run from the `crypto` worktree root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/branch_sage_runtime_info.txt
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/branch_gates.py > research/p192_extension_cover_20261012/branch_gates.log 2>&1
(cd research/p192_extension_cover_20261012/lean && lake build CrossRatio > lean-build.log 2>&1)
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/small_field_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/small_field_pipeline.py > research/p192_extension_cover_20261012/small_field_pipeline.log 2>&1
```

The next useful implementation is a literal \(\mathbb F_{p^2}\) nonhyperelliptic descended cover and a nontrivial \(\mathbb F_{p^3}\) Scholten twist, each with explicit conorm-norm replay. Those constructions must precede any relation solver or full cost. A bounded small-field end-to-end control can then exercise factor-base construction, ordinary relation yield, rank, final linear algebra, and target descent before any P-192 extrapolation.

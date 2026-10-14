# P-192 cover branches: direct obstruction, descended Type I model, and controls

**Result.** Standard NIST P-192 has no direct Joux–Vitse hyperelliptic genus-3 model over the tower \(\mathbb F_{p^6}/\mathbb F_{p^2}/\mathbb F_p\), including the nonsplit-quadratic-root branch. The earlier cross-ratio norm excludes the linear and split-quadratic cases; the nonsplit branch has \(j\notin\mathbb F_p\), while P-192 has \(j\in\mathbb F_p\). The nonhyperelliptic branch is different: we construct an explicit Type I genus-3 cover over \(\mathbb F_p\), hence over \(\mathbb F_{p^2}\), whose degree-two elliptic quotient becomes P-192 over \(\mathbb F_{p^6}\). Its P-192 generator has an explicit quartic coordinate image and two checked cover points. The resulting genus-three class reduces on its \(r\)-primary subgroup to three elliptic quotient coordinates. Its three-step \(p^2\)-Frobenius norm has quotient vector \((-2G,-2G,-2G)\), proving a nonzero \(\mathbb F_{p^2}\)-rational class of exact order \(r\). We also construct a nonfixed-coefficient Scholten genus-2 descent over \(\mathbb F_{p^3}\) and extract a nonzero P-192 \(r\)-torsion class as explicit lower-field Mumford coordinates. A frozen compact-base probe verifies 16 further \(r\)-torsion classes but yields zero three-sum relations in 64 ordinary queries. A separate \(\mathbb F_{23^3}\) instance of the same Scholten descent completes relation collection, full rank, unseen-target recovery, and same-point rho. A scalable P-192 relation method remains the principal implementation gate.

The exact inputs are [in the earlier report](REPORT.md), and the executable receipt is [`branch_gates.py`](branch_gates.py) with [JSON](branch_gates.json), [Sage runtime](branch_sage_runtime_info.txt), and [run log](branch_gates.log). The subgroup is the standard 192-bit prime-order P-192 group; \(p=6277101735386680763835789423207666416083908700390324961279\), \(r=6277101735386680763835789423176059013767194773182842284081\).

The distinction between a subgroup-preserving map and an effective descent is central here. [Momose–Chao, §1](https://eprint.iacr.org/2009/236.pdf) formulates the favorable cover case with an induced isogeny from the cover Jacobian to the Weil restriction of the elliptic curve. The Type I construction follows the paper's descended branch and supplies an explicit affine model and point-level divisor; genus-three Jacobian class arithmetic still needs implementation. The Scholten branch supplies a concrete lower-field subgroup class and a bounded relation probe, enabling a precise test of the first compact-base policy.

## Complete direct hyperelliptic genus-3 exclusion

Let \(p>3\), \(q=p^2\), \(K=\mathbb F_{p^6}\), and let \(d\) be a nonsquare in \(\mathbb F_q\). For a nonsingular quartic in the nonsplit branch
\[
 C_\alpha:y^2=(x^2-d)(x-\alpha)(x-\alpha^q),\qquad \alpha\in K\setminus\mathbb F_q,
\]
its elliptic \(j\)-invariant does not lie in \(\mathbb F_p\). This strengthens the previous P-192 screen from one split-linear-root subfamily to the complete **direct** hyperelliptic genus-3 family of [Joux–Vitse, §4.1 and Appendix A](https://iacr.org/archive/eurocrypt2012/72370010/72370010.pdf). It does not decide whether a different curve in P-192's \(\mathbb F_{p^6}\) isogeny class has such a model.

**Proof.** Choose \(\beta^2=d\) in \(L=\mathbb F_{p^{12}}\). The branch-root cross-ratio gives \(\lambda\in\mu_N\setminus\{1\}\), where \(N=p^4-p^2+1\), and
\(j=256(\lambda^2-\lambda+1)^3/[\lambda^2(\lambda-1)^2]\). These identities and the constructive Cayley parameterization are proved in the existing [two-branch theorem](../iso1_weak_classes_20261007/two_branch_20261009/REPORT.tex). Since \(N\mid p^6+1\), \(\lambda^{p^6}=\lambda^{-1}\). Among the six Legendre transforms with equal \(j\), only \(\lambda\) and \(\lambda^{-1}\) remain in \(\mu_N\): every other transform has norm to \(K\) equal to \(2-(\lambda+\lambda^{-1})\) or its inverse. Norm one would force \(\lambda^2-\lambda+1=0\), hence \(\lambda^6=1\). This is impossible because \(N\) is odd and \(N\equiv1\pmod3\). If \(j\in\mathbb F_{p^d}\), then \(j(\lambda^{p^d})=j(\lambda)\), so \(\lambda^{p^d}=\lambda\) or \(\lambda^{-1}\). For \(d=1,2\), \(\gcd(N,p^d-1)=\gcd(N,p^d+1)=1\); for \(d=2\), the only possible common divisor of \(N\) and \(p^2+1\) is 3, and \(N\equiv1\pmod3\). Hence either equality forces \(\lambda=1\). For \(d=3\), either equality implies \(\lambda^{p^6}=\lambda\), contradicting \(\lambda^{p^6}=\lambda^{-1}\) unless \(\lambda=1\). Thus the degree of \(j\) over \(\mathbb F_p\) is exactly six. The cleared six-factor identity is kernel-checked in Lean as `P192ExtensionCover.legendre_j_difference_factor`. The remaining finite-field torus and norm bridge is written here and in the prior theorem supplement, but is not yet a full Lean theorem.

A linear \(h\), or split squarefree quadratic \(h\), reduces by an \(\mathbb F_q\) coordinate change to the cubic family in [the earlier report](REPORT.md); its required norm-one cross-ratio conflicts with P-192's norm \(-1\). Thus the direct family is exhausted. The exact P-192 \(j\)-invariant is \(6234286251230310114240839169629130138801351179850969208331\) in \(\mathbb F_p\); the run also checked it differs from 0 and 1728.

## Descended Type I nonhyperelliptic genus-3 cover

Let \(F=\mathbb F_{p^3}=\mathbb F_p(\zeta)\), where \(\zeta\) is a root of P-192's 2-division cubic, and put \(q=p^2\). The Legendre parameter \(\lambda=(\zeta^{p^2}-\zeta)/(\zeta^p-\zeta)\) has \(N_{\mathbb F_{p^6}/\mathbb F_{p^2}}(\lambda)=-1\). [Momose–Chao, Lemmas 7–8](https://eprint.iacr.org/2009/236.pdf) reduce their Type I test to the square class of
\[
\Delta_0=(\operatorname{Tr}_{\mathbb F_{p^6}/\mathbb F_{p^2}}\lambda-1)^2-4N(\lambda)
         =(\operatorname{Tr}\lambda-1)^2+4.
\]
Here \(\operatorname{Tr}\lambda\in\mathbb F_p\), so \(\Delta_0\in\mathbb F_p\) and every nonzero such element is a square in \(\mathbb F_{p^2}\): Euler's criterion gives \(a^{(p^2-1)/2}=(a^{p-1})^{(p+1)/2}=1\) for \(a\in\mathbb F_p^\times\). The exact instance gives \(\Delta_0=110791372646815449681899501787236524203341707856604462371\ne0\); it is already a square in \(\mathbb F_p\). The [checked Type I receipt](type_i_gate.json), [Sage runtime](type_i_runtime_info.json), and [SHA-256 list](type_i_receipts.sha256) fix this result. The code solves the paper's Frobenius–Möbius fixed-point quadratic and finds two valid \(\beta\in F\). For the recorded choice, the six elements \(\zeta^{p^i},\beta^{p^i}\), \(i=0,1,2\), are distinct, the defining cross-ratio replays to \(\lambda\), and the quartic
\[
 E_0: Y^2=Q(X)=(X-\zeta)(X-\zeta^q)(X-\beta)(X-\beta^q)
\]
has P-192's \(j\)-invariant. As \(j\ne0,1728\), the two elliptic curves differ over \(F\) by at most a quadratic twist, which becomes trivial over \(\mathbb F_{p^6}\).

The receipt also gives four sparse polynomial equations for a descended affine model. Write \(Y=Y_0+Y_1\zeta+Y_2\zeta^2\), with \(X,Y_i\) over \(\mathbb F_p\). Set the three \(\mathbb F_p\)-coefficients of \(Y^2-Q(X)\) to zero, and add
\[
N_{F/\mathbb F_p}(Y)=N_{F/\mathbb F_p}(X-\zeta)N_{F/\mathbb F_p}(X-\beta).
\]
Take the smooth projective normalization of this affine component. After algebraic closure, the three coefficient equations become the three conjugate quadratic covers. Their product is the square of the six-root polynomial, and the norm equation selects one connected degree-four component. It has six simple branch points, so tame Riemann–Hurwitz gives genus three. Projection to \(E_0\) has degree two after base change to \(F\); after base change to \(\mathbb F_{p^2}\), this is the required \(\mathbb F_{p^6}/\mathbb F_{p^2}\) cover of P-192. Because \(\Delta_0\ne0\), [Momose–Chao, Lemma 8](https://eprint.iacr.org/2009/236.pdf) identifies the smooth genus-3 curve as nonhyperelliptic. This is an explicit descended affine model and a source-backed geometric classification. The P-192 subgroup class is reduced to elliptic quotient coordinates below; a plane-quartic elimination and generic divisor normal form remain separate implementation steps.

Over \(\mathbb F_{p^6}\), the degree-four cover has the Klein-four involutions and three elliptic quotients \(E_0,E_0^q,E_0^{q^2}\). Its quotient by the full group is \(\mathbb P^1\), so the three nontrivial character spaces exhaust \(H^1(C)\). The product of quotient pushforwards is therefore an isogeny \(\operatorname{Jac}(C)\to E_0\times E_0^q\times E_0^{q^2}\), with a 2-power kernel; Galois permutes the factors, giving the Weil-restriction isogeny over \(\mathbb F_{p^2}\). Thus the odd-order P-192 subgroup is preserved. The point-level divisor below gives a concrete conorm and Frobenius norm, and their three quotient coordinates provide the explicit subgroup reduction and norm-order check.

For every \(x\in\mathbb F_p\), put \(U=(x-\zeta)(x-\beta)\). Then \(Q(x)=U^{q+1}\), a square because \(q\) is odd. In fact \(Y=U^{(q+1)/2}\) solves the norm equation too: \(N(U)\in\mathbb F_p\) and \((p^2+1)/2\equiv1\pmod{p-1}\). The other square root has the opposite norm. Thus the policy \(x\in\mathbb F_p\) supplies **exactly \(p\) distinct affine candidate points** on this descended cover. The frozen 256-coordinate sample therefore passes all 256, as checked by the [Sage run](type_i_gate.log) against the same sample digest as the earlier controls. This is a point-count statement for the proposed base policy; the subgroup-usable subset, relation yield, rank, and descent cost have not been measured. Lean checks the algebraic square identity as `type_i_subfield_rhs_square`; the finite-field Frobenius and geometric bridges remain to be formalized.

The first automated \(\beta\)-root attempt stopped at a PARI variable-name error before producing a model; its [raw log](type_i_gate_attempt1.log) is retained. The successful run used the quadratic formula and verified both roots by Frobenius and cross-ratio replay.

### Explicit Type I generator fibers

The [quartic generator map](type_i_generator_map.py) searches the four choices of a quartic root \(\rho\) to send to infinity and the six assignments of the other roots to P-192's three 2-division roots. It finds four valid assignments and fixes one Möbius coordinate
\[
x=C+\frac{D}{X-\rho},\qquad
\kappa^2Q(X)=\bigl(x^3+ax+b\bigr)(X-\rho)^4,
\qquad y=\frac{\kappa Y}{(X-\rho)^2}.
\]
The [receipt](type_i_generator_map.json) retains \(C,D,\rho,\kappa\) and the standard generator's \((X,Y)\in\mathbb F_{p^6}^2\). Sage checks the cleared polynomial identity exactly over \(\mathbb F_{p^3}[X]\), the inverse image and both curve equations, and three further coordinate replays; both additional conjugate quartic right-hand sides are squares at the generator image.

The [divisor transfer program](type_i_divisor_transfer.py) computes the two points \(P_\pm(G)\) and two points \(P_\pm(2G)\) on the four-equation descended model. For each quotient point, it chooses a square root of the first conjugate quartic, obtains the third from the norm equation, and inverts the three-by-three Vandermonde matrix to recover \((Y_0,Y_1,Y_2)\). It verifies all three quadratic equations and the norm equation for each point and every \(p^2\)-Frobenius translate. The explicit degree-zero conorm divisor
\[
P_+(G)+P_-(G)-P_+(2G)-P_-(2G)
\]
represents the pullback of \(G-2G=-G\). The updated certificate checks that each displayed fiber point maps back to its specified P-192 point and that the pushforward of the formal divisor equals \(-2G\ne0\), killed by the odd prime subgroup order \(r\). Since the cover has degree two, \(\pi_*\pi^*=[2]\); pullback is additive, so this conorm class is nonzero and killed by \(r\), hence has exact order \(r\). The group implication and prime-order corollary are kernel-checked as `P192ExtensionCover.conorm_nonzero_torsion_of_push_pull_two` and `P192ExtensionCover.conorm_exact_prime_order_of_push_pull_two`; their geometric push-pull hypothesis follows from the degree-two quotient. Summing the three \(p^2\)-Frobenius translates gives a point-level divisor defined over \(\mathbb F_{p^2}\). The [JSON certificate](type_i_divisor_transfer.json) stores the four initial point coordinates; the [raw log](type_i_divisor_transfer.log) and [Sage runtime](type_i_divisor_runtime_info.json) record the execution.

### Genus-three quotient-vector reduction and Frobenius norm

Let \(\pi_i:C\to E_i\), \(i=0,1,2\), be the three elliptic quotients from the Klein-four action on the Type I cover, and identify each conjugate \(E_i\) with P-192 over \(\mathbb F_{p^6}\) by the \(p^i\)-conjugate of the checked Möbius map. The map \(\Phi=(\pi_{0*},\pi_{1*},\pi_{2*})\) gives an \(r\)-primary class representation: the character identity \(\sum_i\pi_i^*\pi_{i*}=[2]\) on \(\operatorname{Jac}(C)\) makes \(\Phi\) injective on odd-order classes. Indeed \(\pi_i^*\pi_{i*}=1+\tau_i\) for the three nontrivial Klein involutions, while \(1+\tau_0+\tau_1+\tau_2=0\) because the quotient by the full group has genus zero. This is the same Weil-restriction isogeny setting used by [Momose–Chao, §1](https://eprint.iacr.org/2009/236.pdf). It reduces equality, order, and scalar arithmetic on the P-192 subgroup to three ordinary elliptic groups; a generic plane-quartic divisor normal form is not needed for those operations.

The [quotient computation](type_i_norm_quotients.py) maps both points in each fiber through all three conjugate quartic quotients, then repeats after each \(q=p^2\) Frobenius step. The complete frozen vectors are
\[
\Phi(D)=(-2G,0,0),\qquad
\Phi(\sigma D)=(0,0,-2G),\qquad
\Phi(\sigma^2D)=(0,-2G,0),
\]
so \(\Phi(N)=(-2G,-2G,-2G)\) for \(N=D+\sigma D+\sigma^2D\). Each component is nonzero and killed by the prime \(r\). The formal divisor is fixed by \(\sigma\) after the three-step orbit, hence \(N\in\operatorname{Jac}(C)(\mathbb F_{p^2})\). One nonzero quotient component proves \(N\ne0\); Frobenius preserves \(r\)-torsion, so \(N\) has exact order \(r\). The [producer receipt](type_i_norm_quotients.json), [raw log](type_i_norm_quotients.log), and [runtime](type_i_norm_runtime_info.json) retain the three vectors and source hashes. The first attempt used an elliptic curve over \(\mathbb F_p\) for conjugate points over \(\mathbb F_{p^6}\); its [failed log](type_i_norm_quotients_attempt1.log) is preserved. An [independent verifier](verify_type_i_norm_quotients.py) reconstructs the four cover points from the earlier JSON, rechecks all quartic equations, recomputes all nine quotient components and the norm vector, and checks the source/receipt hashes; its [certificate](type_i_norm_verification.json), [log](type_i_norm_verification.log), and [runtime](type_i_norm_verify_runtime_info.json) are retained.

Lean checks the fixedness, nonzero projection, and exact prime order of a three-step norm as `P192ExtensionCover.cyclic_three_norm_exact_prime_order` under the independently checked component facts. It also checks `P192ExtensionCover.three_quotients_faithful_on_odd_torsion`: under the Klein-four push-pull identity, equal quotient vectors determine the same odd-order class. The particular genus-three quotient maps and their Frobenius permutation are Sage-checked, while the geometric derivation of that push-pull identity and finite-field instance remain written mathematical steps. This quotient-vector representation solves the current order and nonvanishing question; constructing a generic plane-quartic divisor reduction and an effective relation collector remain separate implementation tasks.

## Same-field nonhyperelliptic genus-3 control

For \(E:y^2=x^3+ax+b\), define the projective plane quartic
\[
 C: Z^4=X^3W+aXW^3+bW^4,
 \qquad \pi(x,z)=(x,z^2).
\]
It is smooth because, in the affine chart, simultaneous vanishing of its derivatives forces \(z=0\) and a common root of \(f(x)=x^3+ax+b\) and \(f'(x)\); the previous Sage receipt checked \(f\) squarefree. At infinity its sole point \([1:0:0]\) has nonzero \(W\)-derivative. A smooth plane quartic has genus 3 and is nonhyperelliptic. The function-field extension adjoining a square root of \(y\) has degree two. For a degree-two cover, \(\pi_*\pi^*=[2]\), so its pullback is injective on P-192's odd-order \(r\)-subgroup. The polynomial map identity is checked in Lean as `P192ExtensionCover.quartic_projection_identity`.

This construction is defined already over \(\mathbb F_p\), with its degree-two map to P-192 defined over the same field. It establishes a direct subgroup-preserving control distinct from the Type I cover above. A literal plane-quartic Jacobian class for \(G\) and a Diem relation solver remain to be built for this control.

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

This is the fixed-coefficient Scholten specialization: the cover and map are defined over \(\mathbb F_{p^3}\), and P-192 itself is defined over \(\mathbb F_p\). The explicit coefficients in the polynomial basis \(1,\zeta,\zeta^2\) are in the JSON receipt. The route demonstrates literal geometry and subgroup preservation. The next section constructs a nonfixed-coefficient descent; relation arithmetic remains open for both models.

## Nonfixed-coefficient Scholten descent

The [new construction](scholten_twist.py) removes the fixed-coefficient restriction at the level of the cover equation. Retain \(A,B\in\mathbb F_{p^3}\) above, choose a primitive cube root \(\omega\in\mathbb F_{p^6}\), and let \(\sigma:z\mapsto z^{p^3}\). Since \(p^3\equiv2\pmod3\), \(\sigma(\omega)=\omega^2\). The cubic
\[
f_\omega(t)=At^3+B\omega^2t^2+B\omega t+A
\]
has a coefficient \(B\omega^2\) that is not fixed by \(\sigma\), and it satisfies \(f_\omega^\sigma(t)=t^3f_\omega(1/t)\). Its genus-2 double cover over \(\mathbb F_{p^6}\) is \(v^2=f_\omega(u^2)\). The semilinear involution \((u,v)\mapsto(1/u,v/u^3)\) supplies descent data. With \(\eta=\omega-\omega^2\), so \(\eta^2=-3\) and \(\sigma(\eta)=-\eta\), the invariant coordinates
\[
T=\eta\frac{u-1}{u+1},\qquad W=\frac{v}{(u+1)^3}
\]
give the following explicit sextic over \(\mathbb F_{p^3}\):
\[
W^2=\frac{A(\eta+T)^6+B\omega^2(\eta+T)^4(\eta-T)^2
             +B\omega(\eta+T)^2(\eta-T)^4+A(\eta-T)^6}{(2\eta)^6}.
\]
The [Sage receipt](scholten_twist.json), [runtime](scholten_twist_runtime_info.json), and [SHA-256 list](scholten_twist_receipts.sha256) verify every coefficient is fixed by \(\sigma\), the sextic is squarefree and genus two, and the rational map obtained from \(u=(\eta+T)/(\eta-T)\), \(v=W(2\eta/(\eta-T))^3\), and \(x=\zeta+s(\omega u^2+1)\) satisfies P-192's equation. The same frozen 256 \(T\in\mathbb F_p\) values give 126 coordinates with a point over \(\mathbb F_{p^3}\). This supplies a literal lower-field Scholten model with nonfixed coefficient and an explicit map. It may still be isomorphic over the lower field to another specialization; no distinct-isomorphism-class assertion is needed for the cover result. Lean checks the reciprocal cubic identity as `scholten_reciprocal_identity`; Sage checks the fixed coefficients, smoothness, and map, while the literal Jacobian transfer remains a separate gate.

Over \(\mathbb F_{p^6}\), the two involutions \((u,v)\mapsto(-u,v)\) and \((u,v)\mapsto(-u,-v)\) have elliptic quotients \(y^2=f_\omega(t)\) and \(y^2=t^3f_\omega(1/t)=f_\omega^\sigma(t)\). Their pullbacks span the genus-2 Jacobian up to a 2-power kernel. The semilinear descent swaps the quotients, so the product is the expected Weil-restriction isogeny over \(\mathbb F_{p^3}\), preserving odd-order subgroups abstractly.

The [literal transfer run](scholten_transfer.py) checks the standard P-192 generator \(G=(G_x,G_y)\). Put \(t_G=(G_x-\zeta)/s-1\). On the descended sextic, its two preimages have Mumford polynomial
\[
U_G(T)=\omega(\eta+T)^2-t_G(\eta-T)^2,
\qquad V_G(T)=G_y(\eta-T)^3/(2\eta)^3\pmod{U_G},
\]
with \(U_G\) made monic. The code verifies \(V_G^2\equiv h(T)\pmod{U_G}\) and constructs \(D=[U_G,V_G]\) in the genus-2 Jacobian over \(\mathbb F_{p^6}\). Sage checks \(D\ne0\), \([r]D=0\), the same statements for its \(p^3\)-Frobenius conjugate \(D^\sigma\), and \(N=D+D^\sigma\ne0\), \([r]N=0\). As \(\sigma^2=1\), \(N\) is \(\mathbb F_{p^3}\)-rational as a class. The [JSON certificate](scholten_transfer.json), [raw log](scholten_transfer.log), [Sage runtime](scholten_transfer_runtime_info.json), and [SHA-256 list](scholten_transfer_receipts.sha256) retain this replay.

The [lower-field extraction](scholten_lower_field.py) writes every coefficient of the reduced Mumford pair \([U_N,V_N]\) in the basis \(1,\zeta,\zeta^2\) of \(\mathbb F_{p^3}/\mathbb F_p\). It solves a six-by-three prime-field embedding matrix for each coefficient, verifies the base change equals the prior \(N\), constructs \([U_N,V_N]\) afresh in \(\operatorname{Jac}(H)(\mathbb F_{p^3})\), and checks \(N\ne0\), \([r]N=0\). The coefficients and source hash are in the [JSON receipt](scholten_lower_field.json), with [raw log](scholten_lower_field.log), [runtime](scholten_lower_field_runtime_info.json), and preserved [failed accessor](scholten_lower_field_attempt1.log) and [serialization](scholten_lower_field_attempt2.log) runs.

The [bounded split-support probe](scholten_relation_probe.py) takes the exact norm class \(N\) and computes the 256 successive nonzero multiples \(aN\), \(a=1,\ldots,256\). All have reduced Mumford degree two. Exactly 128 of their support polynomials split over \(\mathbb F_{p^3}\); none have both support abscissae in \(\mathbb F_p\). The [JSON](scholten_relation_probe.json), [raw log](scholten_relation_probe.log), [runtime](scholten_relation_runtime_info.json), and [SHA-256 list](scholten_relation_receipts.sha256) fix the query law and counts. This is a finite deterministic split-support diagnostic, not a natural-query rate estimate or a verified relation-matrix row. It demonstrates that the 126/256 eligible \(T\in\mathbb F_p\) cover points do not by themselves yield a subfield-base relation on these 256 queries.

The [compact-base probe](scholten_compact_base.py) then pairs each of 16 further eligible \(T\in\mathbb F_p\) points with one fixed anchor to make valid degree-two Jacobian divisors over \(\mathbb F_{p^3}\). It projects by \(h=\#E(\mathbb F_{p^6})/r\), where \(\#E(\mathbb F_{p^6})\) comes from the exact Frobenius recurrence, and independently checks that every resulting class is distinct, nonzero, and killed by \(r\). The base has 16 subgroup-usable classes and a complete 136-entry pair index. For 64 fixed MT19937 `randrange(1,r)` outputs \(a\), it tests every three-summand equality \(aN=Q_i+Q_j+Q_k\); the frozen run has zero verified relations and rank zero. The [JSON](scholten_compact_base.json), [raw log](scholten_compact_base.log), and [runtime](scholten_compact_base_runtime_info.json) fix the base selection, query digest, source hash, and phase times. The first [failed run](scholten_compact_base_attempt1.log) recorded Sage's even-degree Mumford requirement before the anchor construction was added.

For a fixed base of size \(B\), at most \(\binom{B+2}{3}\) unordered three-sums with repetition can hit distinct nonzero targets. Under an ideal independent uniform nonzero-target law, the per-query hit probability is at most \(\binom{B+2}{3}/(r-1)\), even if every sum is useful. At \(B=16\), this gives \(816/(r-1)\) per query and \(52224/(r-1)<2^{-176}\) for any hit among 64 ideal draws by the union bound. The frozen MT19937 sequence above is deterministic; this probability statement describes the policy under the stated ideal law, not uncertainty in the sealed run. A single-query coverage upper bound of at least one half requires \(B\ge26604808713878332877\), whose complete pair index would have at least \(353907923351028236377112095979348715003\) entries. The [exact integer calculation](compact_base_bound.py) and [receipt](compact_base_bound.json) apply to this unguided exact three-sum policy; structured geometric relation algorithms have a different coverage law.

The first Sage execution reached the coefficient-fixedness check but stopped because the quotient-field element lacks a `.polynomial()` accessor; the [failed log](scholten_twist_attempt1.log) is retained. The verified run extracts the quotient representative with `.lift()` and rechecks the full sextic and elliptic map. An initial Jacobian transfer run in that generic quotient-field representation was capped while Sage searched for a quadratic nonresidue; its exact [base64-encoded interrupted trace](scholten_transfer_attempt1.log.b64) is retained (decoded SHA-256 `9cba89cdbfb2be71483f7242e018467f23632a76ef4afac097698706c062e15a`). The successful replay uses Sage's native \(\mathbb F_{p^6}\) field and the same descended sextic.

## Measured gates and complete-cost boundary

The frozen 256-abscissa sample uses Python MT19937 `randrange(p)` with seed `0x19220261013`; its exact decimal-list SHA-256 is in the JSON receipt. The Sage run completed in 4.487 s on this shared host. These timings are diagnostics, with no isolated CPU timing ratio.

| Branch | Exact direct geometry | Transferred \(r\)-class | 256 abscissae with candidate points | Relation yield / rank / LA / target descent | Full one-target cost |
| --- | --- | --- | ---: | --- | --- |
| Nonsplit quadratic genus 3 | Excluded for P-192 by \(j\)-degree | Not applicable | Not applicable | Not applicable | Not applicable for the direct model |
| Descended Type I nonhyperelliptic genus 3 | Explicit affine normalization over \(\mathbb F_p\), degree-two quotient after base change | Quartic image and both fibers above \(G,2G\) checked; quotient-vector reduction gives norm \((-2G,-2G,-2G)\) of exact order \(r\) over \(\mathbb F_{p^2}\) | 256/256; exactly \(p\) for this \(x\in\mathbb F_p\) policy | Unmeasured | Unknown |
| Same-field nonhyperelliptic genus-3 control | Smooth degree-two plane quartic over \(\mathbb F_p\) | Follows from push-pull; explicit Jacobian arithmetic open | 115; 0.003 s eligibility check | Unmeasured | Unknown |
| Scholten genus 2 | Smooth degree-two cover over \(\mathbb F_{p^3}\) | \(D\ne0, [r]D=0\) checked | 139; 0.005 s eligibility check | Unmeasured | Unknown |
| Nonfixed-coefficient Scholten genus 2 | Explicit sextic over \(\mathbb F_{p^3}\) and rational map to P-192 over \(\mathbb F_{p^6}\) | \(N\) reconstructed as explicit \(\mathbb F_{p^3}\) Mumford class; nonzero and killed by \(r\) | 126/256 eligible; 16 distinct usable projected classes in compact base | 0/64 exact three-sum relations; rank 0/16 for this policy; LA/descent gated | Unknown |

The earlier 115/139 sample tests only square/fourth-power eligibility on the two controls. The Type I model has exactly one point for each \(x\in\mathbb F_p\), while the compact Scholten policy supplies 16 checked subgroup classes and a completed 64-query zero-yield cell. These are different factor-base policies. The original Joux–Vitse policy has \(x\in\mathbb F_p\), so a direct complete scan at P-192 has \(p\) candidate abscissae, about \(2^{192}\). Its reported practical run instead used a roughly 26-bit base prime and a 149-bit subgroup over a degree-six extension [Joux–Vitse, §§3–5](https://iacr.org/archive/eurocrypt2012/72370010/72370010.pdf). The P-192 target-dependent online and same-point rho costs remain unknown under the repository's one-target measurement contract because the compact policy has rank zero and the full-abscissa policy exceeds the available resource envelope.

The exact scan count for this policy is \(p=6277101735386680763835789423207666416083908700390324961279\). Even a one-bit eligibility bitmap indexed by every candidate would require at least
\(\lceil p/8\rceil=784637716923335095479473677900958302010488587548790620160\) bytes, nearly \(2^{189}\) bytes. The candidate-inspection count exceeds \(\sqrt r\) by approximately \(2^{96}\). These are exact work and storage bounds **conditional on that full-abscissa policy**, not measured wall time or a bound on all possible factor-base constructions. The prior 256-coordinate sample estimates only eligibility; it cannot be used to infer relation yield or the number of subgroup-usable base classes.

## End-to-end Jacobian pipeline control

An executable [small-field control](small_field_pipeline.py) exercises ordinary relation collection, useful rank, final modular linear algebra, target descent, and a same-point rho reference on an explicit genus-2 Jacobian. It uses \(E/\mathbb F_{101}:y^2=x^3+x+46\) of prime order 109 and the degree-two cover \(H:v^2=u^6+u^2+46\). Sage verifies \(\#\operatorname{Jac}(H)(\mathbb F_{101})=11663=107\cdot109\), and constructs a nonzero order-109 class \(D=[u^2-2,37]\). Each factor-base class \(P_i\) comes from an affine rational point of \(H\), projected to the order-109 subgroup as \(Q_i=107P_i\). The base contains 32 distinct projected classes. For each ordinary query \(aD\), a pair index searches for an exact three-point Jacobian decomposition \(aD=P_i+P_j+P_k\); multiplying by 107 gives the modular log relation. The code verifies every divisor equality before inserting a row.

| Field and branch | Base classes | Ordinary queries | Verified relations | Independent rows / rank | Fresh-target descent | Same-point rho | Timing status |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| \(\mathbb F_{101}\), same-field genus 2 | 32 | 165 | 63 | 32 / 32 | 7 offset queries; scalar 18 replayed | scalar 18 replayed; 10 walk steps | Shared-host diagnostic |
| \(\mathbb F_{23^3}\), descended Scholten genus 2 | 8 | 12 | 12 | 8 / 8 | 1 query; scalar 15 replayed | same point, scalar 15 replayed; 7 walk steps | Shared-host diagnostic |

The \(\mathbb F_{101}\) run records the deterministic input law, source hash, base digest, each rank-increasing relation witness, phase durations, and target point digest in [JSON](small_field_pipeline.json), with the [Sage runtime receipt](small_field_runtime_info.json), [raw log](small_field_pipeline.log), and [SHA-256 list](small_field_receipts.sha256). Its online target interval was 0.837 s and its matched rho interval 0.000391 s on this shared host.

The new [descended Scholten control](scholten_descended_control.py) uses \(E/\mathbb F_{23}:y^2=x^3+x+4\), an irreducible 2-division cubic, and \(\#E(\mathbb F_{23})=r=29\). It constructs the nonfixed-coefficient genus-two cover over \(\mathbb F_{23^3}\), literally pulls back a generator in \(\mathbb F_{23^6}\), extracts its nonzero \(r\)-torsion norm class \(N\) over \(\mathbb F_{23^3}\), and projects eight point-pair divisors to distinct nonzero \(r\)-torsion base classes. The 12 ordinary queries yield 12 verified relations; eight independent rows give full rank, and final modular linear algebra produces eight base logs checked by direct scalar multiplication. A fresh target is constructed after base-log preparation and recovered in one three-sum query as scalar 15, verified by multiplication. Rho solves that same public point and independently returns 15 in seven walk steps. The final frozen [receipt](scholten_descended_control.json), [raw log](scholten_descended_control.log), and [runtime](scholten_descended_control_runtime_info.json) give the exact cover coefficients, base points, ordinary-query scalars, relation witnesses, target digest, and phase times. A separate [verifier](verify_scholten_descended_control.py) rebuilds the lower-field Jacobian from the receipt, checks the seeded query list, all base projections and relation equalities, re-solves the matrix, and replays the target and rho scalars; its [certificate](scholten_descended_control_verification.json), [log](scholten_descended_control_verification.log), and [runtime](scholten_descended_verify_runtime_info.json) are retained. The first run's [JSON serialization failure](scholten_descended_control_attempt1.log) is preserved.

The frozen descended-control target-dependent intervals are 2.575917 ms for indexed recovery and 10.745583 ms for same-point rho. They exclude setup, base construction, relation collection, and final LA, which are recorded separately in the receipt. CPU isolation was not established, so these timings are algorithmic diagnostics. At P-192 the 16-point compact base has rank zero after its 64-query cap; final LA and a target interval are therefore unavailable for that policy. The control validates the arithmetic pipeline without supplying a P-192 timing extrapolation.

An initial control run failed its scalar replay: target descent omitted multiplication by \(107^{-1}\bmod109\) after cofactor projection. The [failed receipt](small_field_pipeline_attempt1.log) is retained; the corrected formula and the group-projection and scalar identities are checked in Lean as `projected_three_point_relation` and `projected_target_scalar`. This failure illustrates why successful decomposition and full matrix rank alone are not sufficient correctness gates.

## Formal proof map and reproduction

`lean/CrossRatio.lean` builds under Lean 4.34.1 and pinned mathlib revision `d13f23b723b8a846827a245b89c10fc7d3f11612`. Its `#print axioms` output is in `lean/lean-build.log`; the declarations use standard Lean foundations only. The polynomial identity, Scholten coordinate and reciprocal identities, quartic projection and Möbius denominator identities, integer-scalar group projection, field-scalar target formula, odd-Frobenius square identity, exact-prime-order conorm and three-step norm implications, odd-torsion quotient-vector injectivity, fixed/torsion properties of an involution norm, and the exact 816-image three-sum coverage bound are checked. The earlier exact-order failed builds and the failed norm and injectivity attempts are retained as base64-encoded [conorm attempt 1](lean/lean-build-conorm-attempt1.log.b64), [conorm attempt 2](lean/lean-build-conorm-attempt2.log.b64), [norm attempt 1](lean/lean-build-norm-attempt1.log.b64), [injectivity attempt 1](lean/lean-build-faithful-attempt1.log.b64), and [injectivity attempt 2](lean/lean-build-faithful-attempt2.log.b64), preserving their raw bytes without diff whitespace. Sage supplies the P-192 finite-field and quotient-coordinate instances. The Type I finite-field Frobenius instantiation, smooth normalization, Klein-four Jacobian isogeny, exact P-192 elliptic-isomorphism instance, Scholten semilinear descent, Jacobian isogeny, generic plane-quartic divisor reduction, and the probabilistic union-bound specialization are outside those Lean declarations. Full Lean theorems for the direct-family exclusion and descended covers remain open.

| New statement | Lean declaration | Verification boundary |
| --- | --- | --- |
| Project an exact three-point relation by an integer cofactor | `P192ExtensionCover.projected_three_point_relation` | Kernel-checked for any additive commutative group; Jacobian relations separately checked by Sage. |
| Recover the scalar after cofactor projection | `P192ExtensionCover.projected_target_scalar` | Kernel-checked over any field with nonzero cofactor; the modular 109 instance and scalar replay checked by Sage. |
| The Type I quartic RHS is square on fixed subfield coordinates | `P192ExtensionCover.type_i_subfield_rhs_square` | Kernel-checked for a commutative ring given the Frobenius-distributivity equation; the exact P-192 field instance is checked by Sage. |
| The nonfixed Scholten cubic is reciprocal under coefficient conjugation | `P192ExtensionCover.scholten_reciprocal_identity` | Kernel-checked field identity; the P-192 involution and descended sextic are checked by Sage. |
| The Type I Möbius coordinate clears to a quartic polynomial identity | `P192ExtensionCover.mobius_quartic_clear` | Kernel-checked over any field when \(X\ne\rho\); the exact P-192 quartic identity and square twist are checked by Sage. |
| The explicit Type I conorm has exact order \(r\) | `P192ExtensionCover.conorm_nonzero_torsion_of_push_pull_two`, `P192ExtensionCover.conorm_exact_prime_order_of_push_pull_two` | Kernel-checked nonzero, torsion, and exact-prime-order implications under the degree-two push-pull identity; Sage checks the exact P-192 fibers, their quotient images, primality of \(r\), and nonzero \(-2G\). |
| The three-step Type I Frobenius norm is fixed and has exact order \(r\) | `P192ExtensionCover.cyclic_three_norm_exact_prime_order` | Kernel-checked group implication given a 3-cycle, prime torsion, and one nonzero quotient component; Sage and an independent verifier check the concrete component values \((-2G,0,0),(0,0,-2G),(0,-2G,0)\). The finite-field quotient-map instance and Klein-four character identity remain mathematical proof obligations outside Lean. |
| Three elliptic quotient coordinates determine an odd-order genus-three class | `P192ExtensionCover.three_quotients_faithful_on_odd_torsion` | Kernel-checked injectivity given the Klein-four push-pull sum \(\sum_i\pi_i^*\pi_{i*}=[2]\) and an explicit odd-order hypothesis; deriving that sum for the concrete curve remains a written geometric argument. |
| A conorm homomorphism preserves integer torsion | `P192ExtensionCover.conorm_preserves_torsion` | Kernel-checked for any additive commutative groups and additive homomorphism; geometric conorm is supplied by the cover construction. |
| The conjugate-pair norm is fixed by an involution and preserves torsion | `P192ExtensionCover.involution_norm_fixed`, `P192ExtensionCover.involution_norm_torsion` | Kernel-checked group implications; Sage verifies the exact P-192 Mumford class and Frobenius action. |
| A 16-point unguided three-sum base covers at most 816 target classes | `P192ExtensionCover.canonicalTriples16_card`, `P192ExtensionCover.three_sum_support_bound` | Kernel-checked finite triple enumeration and image-cardinality bound; uniform-query union bound and its exact P-192 fraction are checked separately by integer script. |
| Exact P-192 direct-family exclusion and descended-cover isogenies | Open | The finite-field torus, normalization genera, elliptic quotients, and Weil-restriction isogenies are not formalized in Lean. |

Run from the `crypto` worktree root:

The [continuation SHA-256 manifest](continuation_receipts.sha256) seals the revised reports, scripts, outputs, Sage runtime receipts, preserved failed runs, and Lean build log. All JSON source and dependency hashes were checked against their current files before sealing.

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/branch_sage_runtime_info.txt
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/branch_gates.py > research/p192_extension_cover_20261012/branch_gates.log 2>&1
(cd research/p192_extension_cover_20261012/lean && lake build CrossRatio > lean-build.log 2>&1)
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/small_field_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/small_field_pipeline.py > research/p192_extension_cover_20261012/small_field_pipeline.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/type_i_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/type_i_gate.py > research/p192_extension_cover_20261012/type_i_gate.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/scholten_twist_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/scholten_twist.py > research/p192_extension_cover_20261012/scholten_twist.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/scholten_transfer_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/scholten_transfer.py > research/p192_extension_cover_20261012/scholten_transfer.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/scholten_relation_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/scholten_relation_probe.py > research/p192_extension_cover_20261012/scholten_relation_probe.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/scholten_lower_field_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/scholten_lower_field.py > research/p192_extension_cover_20261012/scholten_lower_field.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/scholten_compact_base_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/scholten_compact_base.py > research/p192_extension_cover_20261012/scholten_compact_base.log 2>&1
python3 research/p192_extension_cover_20261012/compact_base_bound.py > research/p192_extension_cover_20261012/compact_base_bound.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/type_i_generator_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/type_i_generator_map.py > research/p192_extension_cover_20261012/type_i_generator_map.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/type_i_divisor_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/type_i_divisor_transfer.py > research/p192_extension_cover_20261012/type_i_divisor_transfer.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/scholten_descended_control_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/scholten_descended_control.py > research/p192_extension_cover_20261012/scholten_descended_control.log 2>&1
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/scholten_descended_verify_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python -u research/p192_extension_cover_20261012/verify_scholten_descended_control.py > research/p192_extension_cover_20261012/scholten_descended_control_verification.log 2>&1
shasum -a 256 -c research/p192_extension_cover_20261012/continuation_receipts.sha256
```

The Type I subgroup class and its three-step Frobenius norm now have explicit quotient-vector reductions and exact order checks. The next Type I gate is a relation mechanism on this \(\mathbb F_{p^2}\)-rational genus-three class, including a usable factor base, verified yield, rank, target descent, and complete cost; a generic plane-quartic divisor normal form is useful if that mechanism requires it. The next Scholten gate is a structured relation mechanism whose base and query costs can compete with the fixed-base three-sum coverage bound. Its ordinary relation yield, useful rank, final linear algebra, and unseen-target descent need bounded P-192 receipts before a complete-cost comparison. The descended \(\mathbb F_{23^3}\) control checks the entire arithmetic pipeline without transferring its costs to P-192.

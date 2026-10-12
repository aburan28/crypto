# ECC2K-130: a degree-7 bridge between the two descending leaf families

Status: **structural deduction and frozen experiment protocol; no new map or PDP run**. The reviewed inputs are the degree-263 kernel/order certificate [#743](https://github.com/aburan28/crypto/pull/743), the exact normalized dual [#750](https://github.com/aburan28/crypto/pull/750), and the equal-useful-size base comparison [#753](https://github.com/aburan28/crypto/pull/753). This note claims no relation, solver gain, or ECC2K-130 logarithm.

## What descent can and cannot change

Let $E/\mathbb F_q:y^2+xy=x^3+1$, $q=2^{131}$, and let $r$ be the 130-bit prime with $\#E(\mathbb F_q)=4r$. The source has $\operatorname{End}(E)=\mathcal O_K=\mathbb Z[\tau]$, where $\tau^2+\tau+2=0$ and $K=\mathbb Q(\sqrt{-7})$. A descending degree-263 map $\phi:E\to L$ has $\operatorname{End}(L)=\mathcal O_{263}=\mathbb Z+263\mathcal O_K$, discriminant $-7\cdot263^2=-484183$, and the *same* endomorphism algebra $K$. The [kernel/order certificate](../../ecc2k130_endo_ring_263_20260925/README.md) checks all 264 source kernel lines in two saved torsion bases: two horizontal lines and 262 descending lines, the latter in two coordinate-squaring cycles of length 131. The [dual certificate](../../ecc2k130_dual_transport_20260925/README.md) checks the normalized full-point identities $\widehat\phi\phi=[263]$ and $\phi\widehat\phi=[263]$ on all eight saved maps, including twist-kernel exceptions.

Since $263\nmid4r$, $\phi$ is a bijection on the rational groups. On the order-$r$ subgroup its inverse is $[263^{-1}\bmod r]\widehat\phi$. For a finite useful leaf base $B_L$, target $T_L=\phi(T)$, and any summand count $m$, put $B_E=\phi^{-1}(B_L)$. The map

$$
(b_1,\ldots,b_m)\longmapsto(\phi(b_1),\ldots,\phi(b_m))
$$

is a bijection between **full-point** decompositions of $T$ over $B_E$ and of $T_L$ over $B_L$. It preserves repetitions, signs, and cofactor-projected group identities. With corresponding factor-base labels, the relation rows over $\mathbb F_r$ and their ranks are identical. This rules out a *support or rank* improvement caused by the conductor change alone. A native leaf base may differ from a chosen source base, but its exact pullback is an equally supported source base. Construction, representation, coordinate-constraint, and orbit-action costs can still differ.

The two Frobenius notions must not be conflated. The geometric 2-Frobenius $\tau:(x,y)\mapsto(x^2,y^2)$ is a cheap degree-2 endomorphism of the source. On a non-horizontal leaf, coordinate squaring goes to a *different* curve. The $q$-Frobenius fixes every point of $L(\mathbb F_q)$ and provides no 131-fold orbit action there. The transported order-$r$ action still exists as a scalar $\lambda$ with $\lambda^{131}=1\bmod r$, but evaluation on native leaf points must be charged. The measured generic affine implementation takes 25,476 field multiplications, 25,604 squarings, and 193 inversions per leaf action; this is an implementation observation, **not a lower bound**. Source $\tau$ takes two squarings in the same certificate. Forward and reverse 263-kernel exceptions are twist-torsion, not rational points of the public order-$r$ group; #750 checks 2,096 nonzero inputs in each direction.

The [equal-useful-size comparison](../../ecc2k130_factor_base_replication_20260925/README.md) found no repeatable descendant-native gain in its separate degree-7 toy holdouts and zero hits for both 16-point, $m=2$ exact degree-263 lines. Under its frozen uniform-target model, expected hits were at most $6.40\times10^{-36}$ per arm; zero is uninformative about an implicit $m\ge3$ PDP. The older [isogeny-class analysis](RESEARCH_ISOGENY_CLASS_SEARCH.md) proves that the Boolean leading form of the *unrestricted direct Semaev* polynomial is curve-coefficient-independent for every $m$. That narrow theorem does not decide affine refutation, rotated base restrictions, full-point encodings, or measured solver cost.

## A new structural pairing, not a new PDP symmetry

The descending order has class number

$$
h(\mathcal O_{263})=h(\mathcal O_K)\,263
\left(1-\frac{(D_K/263)}{263}\right)=263-1=262,
$$

because $h(\mathcal O_K)=1$ and 263 splits in $K$. The 262 leaf isomorphism classes form a class-group torsor. Seven ramifies in $K$ and is coprime to the leaf conductor. [Sutherland, §§2.8–2.10](https://arxiv.org/pdf/1208.5370) therefore gives exactly one horizontal degree-7 isogeny from each leaf. At the class-number-one source, the ramified 7-edge is a self-loop, as the [degree-search note](RESEARCH_ISOGENY_DEGREE_SEARCH.md) records. On a leaf it is **not** a self-loop: the unique invertible ideal of norm 7 is nonprincipal in $\mathcal O_{263}$. A scalar element has square norm and #743 computes 121,046 as the least non-scalar element norm, so there is no norm-7 element. The norm-7 ideal squares to the principal ideal $(7)$; its class has exact order two. This makes the degree-7 map an involution *on leaf isomorphism classes*, not a degree-7 endomorphism of one fixed leaf.

There is also a constructive source-side explanation. The endomorphism
$\alpha=1+2\tau=\sqrt{-7}$ has degree seven, is separable, and satisfies
$\alpha^2=[-7]$. For a descending cyclic kernel $C\subset E[263]$,
$C'=\alpha(C)$ is another descending kernel because $7$ is coprime to
$263$. The universal property of the quotient gives a degree-7 map
$\psi:E/C\to E/C'$ with $\psi\phi_C=\phi_{C'}\alpha$, up to an explicit
codomain isomorphism/orientation. On the *saved twist basis*, coordinate
squaring has matrix $M_{\rm twist}$; source $\tau$ corresponds to
$-M_{\rm twist}$ under the quadratic-twist isomorphism. The exact line
permutation to check is therefore $A=I-2M_{\rm twist}\pmod{263}$, with
$A^2=-7I$. Using $I+2M_{\rm twist}$ would silently apply the wrong
source endomorphism. In the saved seed 20260924 basis, the certificate's
matrix *columns* are $(186,51)$ and $(73,78)$, so the row matrix of
$A$ is $\left(\begin{smallmatrix}155&117\\161&108\end{smallmatrix}\right)$.
Direct modular arithmetic gives $A[1,0]=[1,74]$ as a projective line,
$M_{\rm twist}^{23}[1,4]=[1,74]$, and no other exponent in $0..130$
works. This freezes the predicted conjugacy exponent **23** before any
explicit 7-map is constructed.

The norm-7 class action commutes with coordinate-squaring Galois action. An order-two fixed-point-free permutation commuting with a 131-cycle cannot preserve that odd cycle: within it, any commuting permutation is a power of the cycle, and no nonidentity power has order two. Hence the horizontal 7-edge exchanges the two length-131 leaf families certified by #743. It predicts a testable pairing of the saved nonconjugate representatives [1,0] and [1,4]: the 7-edge codomain from the first is $\mathbb F_q$-isomorphic to a coordinate-squared conjugate of the second for a unique $k\in\{0,\ldots,130\}$. Literal equality of a raw curve coefficient or $j$-invariant is insufficient; full-point $\mathbb F_q$ isomorphism and the twist must be checked. This structure does not restore cheap 2-Frobenius on a fixed leaf and does not create 262 independent challenge trials. The two routes to the paired representation are $E\xrightarrow{\alpha}E\xrightarrow{\phi_{C'}}E/C'$ and $E\xrightarrow{\phi_C}E/C\xrightarrow{\psi}E/C'$. The direct route applies $\alpha$ to both $P,Q$; the bridge route applies $\psi\phi_C$ to both. Their images of unmodified $P,Q$ need not agree without that source action.

There is an exact preflight for constructing the map. The archived order gives $q\equiv4$, $\#E(\mathbb F_q)\equiv2$, and $t=q+1-\#E(\mathbb F_q)\equiv3\pmod7$. Thus the $q$-Frobenius characteristic polynomial on 7-torsion is $X^2-3X+4=(X-5)^2\pmod7$. The source coefficient of $\tau$ in $\tau^{131}$ is $263\cdot146505763881528721$, nonzero modulo 7, so this action is not scalar. Its unique stable ramified-kernel line has eigenvalue 5 of order six; the three signed abscissae should form one degree-three $q$-Frobenius orbit. This is a **derived constraint to verify**, not an observed kernel polynomial. Since $\#L(\mathbb F_q)=4r$ has no factor 7, there is no nonzero rational 7-torsion. The existing rational-abscissa degree-263 Vélu interface cannot simply be reused; an exact cubic kernel polynomial or extension-field construction and independent full-point replay are required.

## Frozen, bounded first gate: certify and price the bridge

This protocol is frozen against origin/main commit 62ef21ec1e083f197593edbe1309c5ddf60b7789, before any bridge construction. Inputs are saved seed 20260924 lines [1,0] and [1,4], the exact field polynomial and public $P,Q,r$ in the challenge fixture, and the reviewed degree-263 forward/dual maps. Current-main input SHA-256 values are:

| Input | SHA-256 |
| --- | --- |
| Saved torsion record | 6f1e7b22f3471214d38ec0d3196c88fe9edaf45791e9c05ea67764831fd75368 |
| Public challenge fixture | 0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f |
| #743 order receipt | bbdff7e901b5d33264e1da95bd0cf1655122d2b1bc42ccce30426a1b24a1640e |
| #750 dual receipt | 4f1d6a5ba85b9c33e07c8af3b42448e10e3633128f3ee9bdffa905fc8cbd2ebc |
| #753 exact replay | 982750f272271a566e48724f470b6bca51d9d3ba159114c7628693b3a2fe9865 |

Any source or input change requires a newly reviewed protocol before its outcome is used for this gate.

1. Compute $A=I-2M_{\rm twist}$ from the frozen torsion certificate and require $A^2=-7I$ and projective equality $A([1,0])=M_{\rm twist}^{k}([1,4])$ for the unique frozen $k=23$ in $[0,130]$. Independently check this line action on a full twist-torsion generator via $\alpha$ under the twist isomorphism. Then construct the unique monic degree-three 7-kernel abscissa polynomial on the [1,0] leaf and verify irreducibility over $\mathbb F_q$, exact point order seven, and Frobenius stability of its cyclic line. Build the full degree-7 Vélu map and normalized dual, descending coefficients from any extension used. Independent arithmetic must validate $\widehat\psi\psi=[7]$ on $P,Q,P+Q,2P$ transported to the first leaf and $\psi\widehat\psi=[7]$ on their images in the codomain. Verify infinity, all six nonzero forward and reverse kernel points in their extension, and off-curve rejection.
2. Compute the codomain model and trace; test every $k=0,\ldots,130$ against coordinate-squared conjugates of the [1,4] leaf. Require that the unique match is the frozen $k=23$ and record an explicit $\mathbb F_q$ isomorphism checked on both public points and eight SHA-256-fixed public coefficient pairs. Use labels leaf-seven-bridge-v1|i|u and leaf-seven-bridge-v1|i|v, with $u=H\bmod r$, $v=1+H\bmod(r-1)$, $i=0,\ldots,7$. Reject a mere $j$-match to the opposite twist. Check $\psi\phi_C(T)=\iota\phi_{C'}\alpha(T)$ on every frozen public control, where $\iota$ is the recorded model isomorphism/orientation, and check $[r]P=[r]Q=O$ after each transport. The two routes' point images need not match without $\alpha$.
3. Compare the **incremental second-leaf stage** after the first degree-263 map is paid: independent second degree-263 setup plus the cheap source $\alpha$ action and transport versus degree-7 kernel discovery, map/dual setup, model normalization, $k$ coordinate squarings, and transport of the same eight targets. Also report cold total including shared torsion discovery and the first degree-263 map. Count all underlying $\mathbb F_q$ multiplications and squarings, including extension arithmetic and inversion internals; report inversion calls separately without double counting. Convert to field-multiplication equivalents with a same-host pre-run calibration of square/multiply over 11 batches of 10,000 SHA-256-fixed nonzero operands (label leaf-seven-cost-v1). Preserve raw counts, CPU, wall and peak RSS. If either calibration has relative median absolute deviation above 5%, or common cold torsion discovery cannot be metered in the same unit, leave the corresponding ratio unset rather than select a favorable conversion.

The structural gate passes only if independent certificates and the unique-isomorphism prediction pass. **Stage-cost promotion** also requires complete incremental cost no more than 0.90 of the direct second-degree-263 route in the frozen multiplication-equivalent unit; otherwise retain the correctness result but stop treating it as a cost lead. Stop after 300 wall seconds or 2 GiB peak RSS per child, preserving a failure, timeout, or OOM receipt. No result from this gate is a PDP, rank, ECDLP, or matched-rho speed claim. A separate reviewed protocol is required before comparing complete implicit $m\ge3$ PDP encodings. That comparison must use a leaf-native base and its exact source pullback on the **same** held-out targets, require per-target witness and rank parity, and charge construction, solver misses, orbit action, transport, and final scalar recovery against same-$Q$ automorphism-aware rho. Only a charged end-to-end win can establish an attack advantage.

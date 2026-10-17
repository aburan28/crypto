# Inversion-free quintic squarefree predicates — derivation and review

Scope: monic quintics over a declared finite field. This is an exact component,
not a new complete solver, no novelty/priority or qualified speed claim.

Write H=X^5+e1 X^4+e2 X^3+e3 X^2+e4 X+e5 in characteristic two.
Put B(T)=T^2+e2 T+e4 and A(T)=e1 T^2+e3 T+e5. Then
H=A(X^2)+X B(X^2), H'=B(X^2). Reduction A mod B is U T+V,
where U=e3+e1 e2 and V=e5+e1 e4. The resultant
D=Res_T(B,U T+V)=V^2+e2 U V+e4 U^2.
Over an algebraic closure every root t of B has a unique square root x,
with multiplicities doubled under T=X^2. At such a root H(x)=A(t).
Consequently Res_X(H',H)=Res_T(B,A)^2=D^2. Thus H is squarefree
iff D != 0, including U=0, V=0, repeated roots of B and every extension degree.
No inversion, root extraction or generic GCD is needed. Seven or fewer field
products/squares suffice after coefficient access. This is a specialization of
classical resultant identities, not established new algebraic theory.

For arbitrary characteristic use a fraction-free polynomial remainder sequence
(H,H'). Pseudo-remainder is ordinary remainder multiplied by nonzero powers
of the divisor's leading coefficient. Over a field this preserves the GCD up
to nonzero scalars. Degree falls each iteration, so the final nonzero polynomial
has degree zero iff H is squarefree. H'=0 correctly returns false for degree5,
including inseparable characteristic-five inputs. Odd extensions are retained.

For squarefree L whose distinct field roots are precisely the allowed domain,
all algebraic roots of H lie in the domain, with repetitions, iff H divides L^5.
In characteristic two H divides L^8 is equivalent. When H squarefree this
reduces to H divides L; otherwise full nilpotence is mandatory. These predicates
classify supplied H; they neither enumerate the hard family nor find its roots.

Independent native oracle: exhaust all monic quintics over F3, F5 and GF4;
reject squarefree iff a monic degree1 or degree2 polynomial has its square
as a divisor. Every repeated irreducible quintic factor has degree <=2, so
this oracle is complete and does not use GCD or the resultant formula.
Support membership is independently compared to literal products over {0,1}.
4,392 squarefree cases and 4,392 support cases pass, including repeated and
inseparable cases. Native Rabin tests cover all9 binary moduli; the odd
quadratic moduli have no base-field root. Sample inverses are checked.
The five NIST characteristic constants are the standard prime formulas,
matching the frozen lab/research/nist_presets.json inputs.

Review: sole-model derivation and independently formulated oracle algorithms;
not an external independent reviewer or Lean-checked theorem. A preparation
compile error left an older executable on disk; its control output is quarantined
as PREPARATION_OLD_BINARY_CONTROLS and is not evidence for the final source.
BUILD_VERIFIED followed by fresh controls succeeded before registration.

Hypothesis: avoid the GCD/inversion cost that erased the previous shortcut.
Prediction: native specialized characteristic-two resultants cost less than
native Euclidean classification on wide binary fields. Generic PRS may regress;
retain all fields and reject blanket adoption on a regression. Two seeds fixed
before comparison; all three methods share native arithmetic. Timings are
three alternating internal rounds, not three cold solver pipelines or F4/F5/F6.
The updated8x goal is separate from, and does not weaken, performance/GATE.md.
Primary-source novelty review is unfinished; PRS/resultants are classical.

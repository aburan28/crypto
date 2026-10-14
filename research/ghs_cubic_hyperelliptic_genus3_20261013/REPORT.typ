#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set text(size: 11pt)

= The hyperelliptic genus-3 GHS cover ($m = 3, r = 5$): explicit $FF_q$ model, Frobenius, and the transfer executed with Cantor arithmetic

2026-10-13 · Conductor T-2463 · `research/ghs_cubic_hyperelliptic_genus3_20261013/` · recorded run `20261014T0029Z` (base commit `12331a9a7`, source sha256 `373c1663…`). *Structural stage diagnostic on toy fields.* No index calculus, no relation collection, no discrete logarithm, no cost or speedup claim. Protocol frozen before the recorded run (`PROTOCOL.md`, one dated amendment written after development runs).

== Question
The cover census (PR \#1759) found a hyperelliptic genus-3 type ($m=3$, $r=5$) of GHS cover for $E\/FF_(q^3)$. The $m=2$ genus-3 covers built in PR \#1761 are plane quartics without class arithmetic in the repository. Can the $m=3, r=5$ curve be written over $FF_q$, does its Frobenius match, and does the transfer $E_0(FF_(q^3)) -> "Jac"(C)(FF_q)$ run as a group map?

== Construction (derived, then verified)
$K = FF_(q^3)$, $p in K without FF_q$, $r_i = p^(q^i)$, $Q in FF_q [x]$ monic squarefree (split $(x-a)(x-b)$ or irreducible), $F_0 = Q(x-r_0)(x-r_1)$, $F_i = F_0^(sigma^i)$, $N = product (x - r_i)$. $C = K(x, y_0, y_1, y_2)$, $y_i^2 = F_i$, descent $y_0 |-> y_1 |-> y_2 |-> y_0$. The quotient $D = C\/iota$ is the conic $(r_2-r_0)X_0^2 + (r_0-r_1)X_1^2 + (r_1-r_2)X_2^2 = 0$ in $X_0 = u_(12)u_(20)$, … with $u_(01) = y_0 y_1\/(Q(x-r_1))$; in trace coordinates it is defined over $FF_q$. A line parametrization $z$ through an $FF_q$-point gives $x = x_n (z) \/ x_d (z)$ and $w = y_0 y_1 y_2 x_d (z) \/ (Q N)$ with $w^2 = h(z) = x_n^2 + c_1 x_n x_d + c_0 x_d^2$ (since $F_0 F_1 F_2 = Q^3 N^2$). Split arm: the unique $FF_q$-point of $D$ over $a$ at $z = infinity$ gives $deg h = 7$. Transfer $phi(P) = sum_(j=0)^2 "Fr"^j (pi_0^* [P] - 4 infinity)$, base $O = (r_0, 0)$. Derived: $pi_(0*) compose phi = [4]$, so $ker phi subset E_0(K)[4]$; the unramified subcover $K(E_0)(u_(12))$ puts $T_1 = [(r_1,0)] - [O]$ in $ker phi$.

== Result (20 instances, $q = 5, 7, 11, 13$; 0 failures)
#table(columns: (auto, 1fr), stroke: 0.4pt, inset: 4pt,
  [*P1 descent*], [$FF_q$ coefficients 20/20; on 200 random points per instance the map lands on $w^2 = h(z)$, recovers $x$, commutes with Frobenius and inverts (4000/4000 each)],
  [*P2 genus 3*], [$deg x(z) = 4$, $h$ squarefree of degree 7 (split) / 8 (conjugate), 20/20],
  [*P3 Frobenius*], [$N_1, N_2, N_3 = q+1, q^2+1, q^3+1-3t_E$ on the $FF_q$ model, equal to the fibre count on the $K$-model; $L(T) = 1 - t_E T^3 + q^3 T^6$, $L(1) = \#E_0(K)$, 20/20; Mumford-pair enumeration gives the same order, 12/12 split],
  [*P4 transfer* (split)], [$phi(P)$ reduced over $FF_q$ for all points; homomorphism 400/400 pairs and 100/100 doublings per instance; $phi(G) != 0$, $ell phi(G) = 0$ (12/12); $|ker phi| in {2, 4}$, $subset E_0(K)[4]$, $T_1 in ker$ (12/12); $pi_(0*) phi(P) = 4P$ for every $P in E_0(K)$ (12/12)],
  [*P5 census*], [$q = 5, 7$: census class of $E_0$ and both conjugates has min genus 3, type $m=3, r=5$ hyperelliptic, same order (11 instances); $q = 11, 13$ not in the census sample. All 167,472 enumerated candidates have $4 | \#E_0(K)$],
)

#figure(image("reports/run-20261014T0029Z/checks.svg", width: 100%), caption: [Run `20261014T0029Z`. Left: $\#"Jac"(C)(FF_q)$ from $L(1)$ (20 instances) and from Mumford enumeration (12 split) against $\#E_0(K)$. Right: transfer checks per split instance; kernel order annotated.])

== Example over $FF_7$
$K = FF_7[t]\/(t^3+3)$, $p = t$, $Q = x(x-1)$, $\#E_0(K) = 328 = 8 dot 41$, $t_E = 16$:
$C: w^2 = 6z^7 + z^6 + 3z^5 + 2z^4 + 2z^3 + 5z^2 + 6z + 6$, $x(z) = (z^3+5z^2+6z+6)\/(z^4+2z^3+5z^2+6z+5)$. $N_k = 8, 50, 296$; $\#"Jac"(C)(FF_7) = 328$. $phi(G) = (z^3+2z^2+5z+4, 4z^2+2z+3)$ with $41 phi(G) = 0$; $ker phi$ cyclic of order 4, image of size 82.

== Lean
`lean/HyperGenus3.lean` (v4.34.1, no Mathlib, no `sorry`): branch-set combinatorics ($r=5$, $m=3$, genus 3, rank-2 subspace with union of size 3, $tau B Delta tau^2 B subset B$, triple combination $= {a,b}$), the conic identities, $F_0F_1F_2 = Q^3N^2$, and $w^2 = h$ multiplied out. Not formalized: function-field descent, Frobenius, L-polynomial, the transfer and $pi_(0*) phi = [4]$ (numerical only).

== Limits
Closes the plane-quartic obstruction *by a different route for this type only*; plane-quartic class arithmetic is still missing. $phi$ is not an isomorphism (index 2 or 4). No transfer on the conjugate arm (degree-8 model, no real-model Cantor). Toy fields, first instances in a fixed order. No cost statement.

# Contract and branch proof

Binary component contract: decide whether every algebraic root, with its
multiplicity, of a supplied monic quintic H lies in the affine domain
D=a+span_F2(b_0,...,b_17). This is a predicate, not root enumeration or an
original norm-family search. No novelty claim. All arithmetic is exact.

Write H=X B(X^2)+A(X^2), B=T^2+e2*T+e4 and
A=e1*T^2+e3*T+e5. Put U=e3+e1*e2, V=e5+e1*e4 and
Delta=V^2+e2*U*V+e4*U^2. Euclidean elimination gives
A+e1*B=U*T+V. Therefore H is squarefree iff Delta!=0.
If Delta=0 and U!=0, t=V/U is a root of B and of A. In a finite
characteristic-two field it has a unique square root x, and (X+x)^2=X^2+t
divides H. Every root of H lies in D iff x belongs to D and every root of
the exact cubic quotient Q=H/(X^2+t) lies in D. Multiplicities, including
triple/quintuple roots or another repeated root in Q, do not change this iff.
If U=0, Delta=0 forces V=0; the implementation takes the exact generic
fallback and does not divide by zero or discard that branch.

Let L_0(Z)=Z, beta_i=L_i(b_i), and
L_(i+1)=L_i^2+beta_i*L_i. Its distinct roots are the prefix subspace,
and beta_i!=0 checks independence. D is the root set of L_d(Z)+L_d(a).
For a monic degree-m polynomial, all roots belong to D iff the remainder
of (L_d+L_d(a))^(2^ceil(log2(m))) modulo H is zero. Each irreducible
factor must be linear with root in D; the exponent covers every multiplicity.
For squarefree H the unsquared remainder alone is necessary and sufficient.
The comparator already uses the same squarefree shortcut.

For the repeated root, x in D iff t=x^2 belongs to D^2. Squaring all
coefficients of L_i gives the annihilator of the squared basis.
For U!=0 maintain P_i=U^(2^i)*L_i^(coefficient square)(V/U), W_i=U^(2^i).
Initially P_0=V, W_0=U. Then
P_(i+1)=P_i^2+beta_i^2*W_i*P_i and W_(i+1)=W_i^2.
The exact affine membership condition is
P_d+L_d(a)^2*W_d=0. This uses no inverse or square root. A false condition
rejects all-support, but DOES NOT assert that H has no roots in D.
A true condition permits an ordinary field inverse to produce t and Q.
The exceptional U=V=0 branch and every squarefree case are preserved.

For odd characteristic p the separate p-linear recurrence is
L_(i+1)=L_i^p-beta_i^(p-1)*L_i.
For rational V/U in a+V_i, set P_i=U^(p^i)*L_i(V/U), W_i=U^(p^i):
P_(i+1)=P_i^p-beta_i^(p-1)*W_i^(p-1)*P_i,
W_(i+1)=W_i^p. Membership iff P_d-L_d(a)*W_d=0.
The native implementation tests this against literal domain products in
small prime fields and odd quadratic extensions. This is not an odd quintic
factorization algorithm. Prime-field Fp has only dimension0/1 over Fp;
this recurrence supplies no proper large prime-field interval domain.

Work: degree<=5 quotient operations and O(d) field operations per branch;
small fixed-degree storage and O(d) domain constants. Input generation uses
outside=a+scale*X^18 with nonzero scale in degree>18 fields, guaranteeing
nonmembership and removing the earlier prefix-probe loop. Complete paid
predicate includes domain setup once per batch, discriminator and quotient
work. Timings separately expose amortized per-domain batches; they are not
cold original-solver measurements. No parameter-family enumeration is removed.

Review: one-model derivation with an independently formulated Hasse-derivative
oracle; not an independent-human or Lean review. Tiny exhaustive GF4 quintics,
sampled GF8/GF16 quintics, outside/repeated/nonrational cases retain all branches.
No unchecked Lean certificate is supplied. Classical finite-field factoring
and linearized subspace polynomials are prior art; priority is not established.

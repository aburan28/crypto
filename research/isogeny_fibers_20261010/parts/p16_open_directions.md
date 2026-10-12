---

## 16. Open directions (2026-10-11)

§15.4 leaves no route to an exponent change, but it does leave measurable
questions where the constants, the feasibility frontier or the defender's
cost could move. Ranked by expected value. Each item names the first
experiment and what its result would mean.

### 16.1 Where kernel saturation pays: the feasibility crossover in m

The affine law (Prop. 15.2) charges a kernel-saturated base |M|^{m−2} in
yield. Symmetrised summation polynomials buy degree reductions that also
grow with m: S₄ drops from degrees [4,4,4,4] to [2,2,2,2,1]. The previous
round found that at n = 17, m = 3, the symmetrised system was the only one
the engine could solve at all (§7.A1-R and the 2026-09-27 contract). So
there is a regime where the penalty buys feasibility.
*Experiment.* For n ∈ {17, 19, 23} and m ∈ {3, 4}, run plain x-frame and
u-frame symmetrised arms under one wall cap. Record solved or unsolved,
S and hit rate. Find the smallest (n, m) where the plain arm is infeasible
and the symmetrised arm is not, and compare its S with the best feasible
m = 2 arm on the same curve.
*Meaning.* If the crossover arm ever beats every feasible arm at the
same n, kernel saturation earns a real role. The law predicts it will
not, because m = 2 arms are 50× cheaper at n = 17. A measured exception
would be the first counterexample to §15.4's row on algebraic oracles.

### 16.2 A fast 4E test as a subgroup check (defender side)

Prop. 14.2 decides P ∈ 4E on y² + xy = x³ + 1 with two traces, one
half-trace and one multiplication. When h = 4, as for ECC2K-130 and
sect233k1, 4E is exactly the prime-order subgroup. That replaces the
[r]P membership check, about n doublings and additions, with O(n)
squarings.
*Experiment.* Implement the test in the binary-field fast path. Validate
it against [r]P on 10⁶ random points for n ∈ {131, 233}, and benchmark it.
*Meaning.* A practical constant-factor win for point validation. The
trace test for halvability is classical; the two-level closed form and its
use as a full subgroup check for cofactor-4 Koblitz curves should be
checked against the literature before any claim.

### 16.3 The vertical security margin

Thm. 15.1 and Prop. 3.2 leave one region open: a conductor gap with a large
prime factor ℓ, where the transfer isogeny exists but costs Õ(√ℓ) to
evaluate after Õ(ℓ²) to find.
*Experiment.* Time kernel finding (ψ_ℓ factor or modular-polynomial root)
and √élu evaluation for ℓ up to 2^20 with the repository's isogeny code.
Fit the exponents, and report for a given class the ℓ above which a
vertical transfer costs more than rho on the target.
*Meaning.* Turns "polynomial in ℓ" into a number per curve: how large the
conductor's prime factors must be before the vertical neighbours are
cryptographically independent.

### 16.4 Algebraic oracles over F_{p²}: the free automorphism fold

On a j = 0 curve over F_{p²} with p ≡ 1 mod 3, the Gaudry base {x ∈ F_p}
is automatically ζ-closed, because ω ∈ F_p. The affine law then predicts
a free ×3 cut in unknowns at unchanged Weil-descent cost. A 3-torsion
saturated alternative is predicted to lose 3× at m = 3 despite its
smaller symmetric system.
*Experiment.* Run both arms through the repository's F_{p²} pipeline
(`glv_invariant_base.rs`, `gls_fp2.rs`) at matched p and wall cap, and check
whether the ζ-fold is already applied.
*Meaning.* Confirms Prop. 15.2 for an algebraic oracle in odd
characteristic, where so far it is only measured for the mitm oracle.

### 16.5 Class-wide multi-target rho, measured

§14.5 predicts that L keys on L isogenous curves cost about √(Lr) in
total once transported.
*Experiment.* Take a toy class with explicit 2- and 3-isogeny paths
(`csd/volcano.py` builds them) and L ∈ {1, 4, 16, 64}. Transport the
targets and run one batch rho through ecbench against L independent runs.
*Meaning.* Turns a folklore argument into a measured guidance line for
protocols that issue per-user curves.

### 16.6 Genus 2: does a larger automorphism group change the law?

Jacobians of y² = x⁵ + 1 carry an order-10 automorphism group and admit
rational (2,2)-isogeny kernels. The law predicts a gain of |A| and a loss
of |M|^{m−2}, nothing more.
*Experiment.* Port the exact census of §15.2 to a genus-2 Jacobian over a
small prime field, with Mumford-coordinate bases.
*Meaning.* Either the law is dimension-free, or higher genus has a
saturation regime the elliptic case lacks.

### 16.7 Theory: isogeny oracles in the algebraic group model

Thm. 15.1 is generic. Index calculus is not, because it reads coordinates.
The natural next statement adds isogeny oracles to an algebraic-group-model
reduction and asks whether the summation-polynomial structure can be
expressed there.
*Goal.* Either a reduction showing that coordinate access plus isogeny
oracles adds nothing beyond coordinate access alone, or an explicit
construction where it does.

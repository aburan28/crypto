# Exact fixed-U affine cell intersection

Contract: nonzero fixed U in K=GF(2^n), independent b_0..b_(d-1),
affine D=a+span_F2(b), and V(z)=v0+sum z_j*v_j, z_j in F2.
Return a compact affine basis of ALL parameter masks z for which V(z)/U
belongs to D^2. This is a necessary repeated-root support condition, not
all original roots, not all-support on its own, and not a bound to the
actual norm-family parameterization. No novelty claim; linear intersection
is classical. U=0 is explicitly outside this component contract.

Condition is V(z)+U*a^2 in W=span_F2(U*b_i^2). Squaring and multiplication
by U are F2-linear bijections, so W has dimension d. Construct a Gaussian
basis of W and a linear quotient projection pi:K->K/W. Solve
sum z_j*pi(v_j)=pi(v0+U*a^2). Column elimination tracks each pivot's
parameter mask. Each dependent column provides a kernel vector whose
highest new parameter bit is unique; hence these vectors are independent.
They span the kernel by rank-nullity. Reducing the RHS either witnesses
inconsistency or produces a particular mask. Its translate by the kernel
is therefore the entire solution set. Returning an affine basis is explicitly
a different output representation from enumerating original root tuples.

If the V map is injective, the surviving parameter dimension is at most d:
the image of the solution differences lies in W. The worker uses the first
n-10 polynomial basis vectors, hence injectivity is guaranteed. d=18 fixes
an upper bound18 for all six fields; generic dimension8 is a prediction,
not a theorem about these seeds or about the original2^73 norm family.
Gaussian work is polynomial in n,k,d with bitset XOR; storage O(n*(n+k))
bits plus the original columns. No 2^k parameter scan is used in the worker.
Tiny independent controls enumerate150 cells and9600 parameter memberships.
Certificate checks use a separately formulated annihilator recurrence.

The affine offset must be squared before evaluating the squared-basis
annihilator. A pre-freeze review corrected a helper that initially squared
the wrong annihilator evaluation; no wide frozen experiment ran that helper.
The final tiny controls test both the Gaussian solution and this verifier
against literal domain products. Old preparation files remain preserved.

Critical open transfer: U may vary nonlinearly in actual hard-family
coefficients. Partitioning by every field value of U would reintroduce an
exponential search and is NOT proposed. Next inspect the actual source map
for naturally fixed-U cells, or derive a bilinear incidence solver retaining
variable U. Retain Delta=0, squarefree, exceptional U=V=0 and source-fiber
coverage separately. This component does not solve the global search gap.

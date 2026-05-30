\\ Check the natural degree-6 model y^2 = prod_{i=1..6}(x - r_i)
\\ where {r_1,...,r_6} = E1 2-torsion UNION E2 2-torsion.
\\
\\ For p=13, b=1, d=2:
\\   E1: y^2=x^3+1 mod 13, 2-torsion at x = 12, 4, 10 (roots of x^3=-1=12)
\\   E2: y^2=x^3+8 mod 13, 2-torsion at x = 7, 8, 11 (roots of x^3=-8=5)
\\
\\ The natural genus-2 model: C: y^2 = (x-12)(x-4)(x-10)(x-7)(x-8)(x-11)
\\ If Howe conditions hold, Jac(C) should be (2,2)-isogenous to E1 x E2.
\\
\\ Run: gp -q chlrs_sextic.gp

default(parisize, 32000000);
p = 13;

E1 = ellinit([0, 1], p); t1 = p + 1 - ellcard(E1);
E2 = ellinit([0, 8], p); t2 = p + 1 - ellcard(E2);
print("t1 = ", t1, "  t2 = ", t2);
print("Expected char poly of E1 x E2: x^4 + ", 2*p - t1^2, "*x^2 + ", p^2);

\\ E1 2-torsion: roots of x^3 + 1 mod 13
\\ E2 2-torsion: roots of x^3 + 8 mod 13

\\ Natural sextic: use Fp coefficients (Mod(1,p) * ...)
h6 = Mod(1,p) * (x-12)*(x-4)*(x-10)*(x-7)*(x-8)*(x-11);
print("Sextic h6(x) = ", h6);
cp6 = hyperellcharpoly(h6);
nj6 = subst(cp6, variable(cp6), 1);
target = (p+1-t1)*(p+1+t1);
print("Char poly of Jac(C) (sextic model): ", cp6);
print("#Jac(C) = ", nj6, "  target = ", target, "  match = ", nj6 == target);

\\ Now try each of the 3 possible Rosenhain models (different partitions into triples)
\\ Partition 1: map E1 2-torsion triple to {0, inf, ...} by sending
\\   12->0, 4->inf (so one E1 point goes to inf, giving degree-5 model)
\\ Apply Mobius T: x -> (x-12)/(x-4) * constant

\\ Alternative: apply an affine transform x -> a*x + b to move one branch
\\ point to 0 and another to 1, then the sextic becomes pentadic.

\\ Actually, PARI can handle degree-6: just verify sextic gives correct Jac.

print("");
print("Try alternative orderings of branch points:");

\\ Rosenhain form: move one branch point to infinity.
\\ Send r_1 -> inf via x -> 1/(x - r_1). The new h becomes:
\\ y^2 = (...) degree 5.

\\ Manual: send 12 -> inf by substituting x -> 1/(x-12)+12 (or just use Mobius).
\\ Simpler: for the degree-5 form, the 5 remaining branch points are 4, 10, 7, 8, 11.
\\ After Mobius T(z) = (z-12)/(z-4) (sending 12->0, 4->inf, need another point->1):
\\   T(10) = (10-12)/(10-4) = -2/6 = -1/3 = -1*3^{-1} mod 13.
\\   3^{-1} mod 13: 3*9=27=1 mod 13, so 3^{-1}=9. T(10) = -9 = 4 mod 13.
\\   T(7) = (7-12)/(7-4) = -5/3 = -5*9 = -45 = -45+4*13 = -45+52 = 7 mod 13.
\\   T(8) = (8-12)/(8-4) = -4/4 = -1 = 12 mod 13.
\\   T(11) = (11-12)/(11-4) = -1/7. 7^{-1} mod 13: 7*2=14=1, so 7^{-1}=2. T(11)=-2=11.
\\ Branch pts after T: {0, inf, 4, 7, 12, 11} - wait, T(10)=4, T(7)=7, T(8)=12, T(11)=11.
\\ We need to scale so that ONE of the images equals 1.
\\ Scale: divide by T(7)=7, so T' = T/7. T'(10)=4/7=4*2=8, T'(7)=1, T'(8)=12/7=12*2=24=11, T'(11)=11/7=11*2=22=9.
\\ Under T'(z)=(z-12)/(7(z-4)):
\\   T'(12)=0, T'(4)=inf, T'(10)=8, T'(7)=1, T'(8)=11, T'(11)=9.
\\ Rosenhain form: y^2 = x*(x-1)*(x-8)*(x-11)*(x-9) (degree-5, branch pts 0,1,inf,8,11,9).
h5b = Mod(1,p)*x*(x-1)*(x-8)*(x-11)*(x-9);
cp5b = hyperellcharpoly(h5b);
nj5b = subst(cp5b, variable(cp5b), 1);
print("Rosenhain (12->0,4->inf,7->1): char poly = ", cp5b, "  #Jac = ", nj5b, "  match = ", nj5b == target);

\\ Also try Rosenhain with E2-pt at inf
\\ Send 7->inf via T(z)=(z-12)/(z-7) (12->0, 7->inf).
\\   T(4)=(4-12)/(4-7)=-8/-3=8/3=8*9=72=7*13+1... 72=5*13+7. T(4)=7*9 mod 13? Let me recompute.
\\   3^{-1}=9 mod 13. T(4)=(-8)*(-3)^{-1}=-8*(-3)^{-1}. (-3)^{-1}=-(3^{-1})=-9=4. T(4)=-8*4=-32=-32+3*13=-32+39=7.
\\   Normalise by T(4)=7: T'=(z-12)/(7(z-7)). T'(4)=1.
\\   T'(10)=(10-12)/(7*(10-7))=-2/(7*3)=-2/21. 21 mod 13=8. -2*8^{-1}=-2*5=-10=3. T'(10)=3.
\\   T'(8)=(8-12)/(7*(8-7))=-4/7=-4*2=-8=5. T'(8)=5.
\\   T'(11)=(11-12)/(7*(11-7))=-1/(7*4)=-1/28. 28 mod 13=2. -1*2^{-1}=-7=6. T'(11)=6.
\\ Rosenhain: 0,1,inf,3,5,6. h = x(x-1)(x-3)(x-5)(x-6).
h5c = Mod(1,p)*x*(x-1)*(x-3)*(x-5)*(x-6);
cp5c = hyperellcharpoly(h5c);
nj5c = subst(cp5c, variable(cp5c), 1);
print("Rosenhain (12->0,4->1,7->inf): char poly = ", cp5c, "  #Jac = ", nj5c, "  match = ", nj5c == target);

print("");
print("Target: x^4 + 22*x^2 + 169 (char poly of E1 x E2)");
print("If none match: the Howe-gluing formula for j=0 curves needs Weil-pairing correction.");

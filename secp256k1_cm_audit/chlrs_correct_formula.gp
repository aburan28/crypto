\\ Find the algebraic relationship between (b, b2, d) and the correct c values.
\\ Context: for p=1009, E1: y^2=x^3+11, E2 (quad. twist): y^2=x^3+515 (d=11, b2=d^3*b)
\\ Correct c: y^2=(x^3+11)(x^3+c) has Weil poly x^4+169*x^2+p^2 for c in {140, 453, 919}.

default(parisize, 256000000); default(timer, 0);

p = 1009; b = 11; d_ns = 11;
b2 = lift(Mod(d_ns^3 * b, p));   \\ = 515
om = lift(Mod(1,p));
{
  my(r); r = polrootsmod(Pol([1,1,1]), p);
  if(#r == 0, error("no cube root of unity"));
  om = lift(r[1])
};
print("p=", p, " b=", b, " b2=", b2, " d=", d_ns, " om=", om);
print("Correct c values: {140, 453, 919}");
print("");

\\ Test candidate formulas for c:
c_vals = [140, 453, 919];
print("---- Candidate formulas for c ----");

\\ 1. c = b^2 / b2 mod p?
f1 = lift(Mod(b^2, p) * Mod(b2, p)^(-1));
print("b^2/b2 mod p = ", f1, "  in {140,453,919}? ", vecsearch(c_vals, f1) > 0);

\\ 2. c = b / d^3 mod p?
f2 = lift(Mod(b, p) * Mod(d_ns^3, p)^(-1));
print("b/d^3 mod p = ", f2, "  in {140,453,919}? ", vecsearch(c_vals, f2) > 0);

\\ 3. c = b2^2 / b mod p?
f3 = lift(Mod(b2^2, p) * Mod(b, p)^(-1));
print("b2^2/b mod p = ", f3, "  in {140,453,919}? ", vecsearch(c_vals, f3) > 0);

\\ 4. c = sqrt(b * b2) mod p?
bb2 = lift(Mod(b*b2, p));
f4 = lift(Mod(bb2, p)^((p+1)/4));
print("sqrt(b*b2) mod p = ", f4, "  (b*b2=", bb2, ")  sq check: ", lift(Mod(f4^2,p)) == bb2);

\\ 5. c = b * om mod p and variants?
f5a = lift(Mod(b * om, p));
f5b = lift(Mod(b * om^2, p));
print("b*om mod p = ", f5a, "  in {140,453,919}? ", vecsearch(c_vals, f5a) > 0);
print("b*om^2 mod p = ", f5b, "  in {140,453,919}? ", vecsearch(c_vals, f5b) > 0);

\\ 6. c = b2 * om or b2 * om^2?
f6a = lift(Mod(b2 * om, p));
f6b = lift(Mod(b2 * om^2, p));
print("b2*om mod p = ", f6a, "  in {140,453,919}? ", vecsearch(c_vals, f6a) > 0);
print("b2*om^2 mod p = ", f6b, "  in {140,453,919}? ", vecsearch(c_vals, f6b) > 0);

\\ 7. c = (b * b2)^(1/3) or its twists?
bb2 = lift(Mod(b * b2, p));
print("");
print("b*b2 mod p = ", bb2);
print("Is b*b2 a cube? ", lift(Mod(bb2,p)^((p-1)/3)) == 1);
if(lift(Mod(bb2,p)^((p-1)/3)) == 1,
  f7 = lift(Mod(bb2,p)^((2*(p-1)/3+1)/3));
  print("(b*b2)^(1/3) = ", f7, "  in {140,453,919}? ", vecsearch(c_vals, f7) > 0)
);

\\ 8. Try: c = b * g^k for various k, where g is primitive root
g = lift(znprimroot(p));
print("");
print("Primitive root g = ", g);
{
  for(k = 0, p-2,
    ck = lift(Mod(b * g^k, p));
    if(vecsearch(c_vals, ck) > 0,
      print("  c = b * g^", k, " = ", ck, "  [g^k = ", lift(Mod(g^k,p)), "]");
      print("    g^k = ", lift(Mod(g^k, p)), " = ?")
    )
  )
}

\\ 9. What are the relationships among the three c values?
print("");
print("---- Relationships among {140, 453, 919} ----");
r01 = lift(Mod(c_vals[2], p) * Mod(c_vals[1], p)^(-1));
r02 = lift(Mod(c_vals[3], p) * Mod(c_vals[1], p)^(-1));
print("453/140 mod p = ", r01, "   vs. om=", om, "  om^2=", lift(Mod(om^2,p)));
print("919/140 mod p = ", r02, "   vs. om=", om, "  om^2=", lift(Mod(om^2,p)));
print("");

\\ 10. Trace formula: Richelot image of y^2=(x^3+a)(x^3+b)
\\ For the (2,2)-isogeny from Jac(C) to a product, the dual curve has
\\ a specific form. For biquadratic genus-2, the Richelot formula gives:
\\ If C: y^2 = (x^3+a)(x^3+b), decompose as y^2 = G1*G2*G3 for quadratics.
\\ But for biquadratic f = x^6 + (a+b)x^3 + ab, the (2,2)-kernel
\\ uses the factorization into quadratics.
\\
\\ For the specific case f = x^6 + S*x^3 + P (where S=a+b, P=a*b),
\\ the discriminant of x^3 + S/2*x^3/2 + P... actually let me try:
\\ x^6 + S*x^3 + P = (x^2 + px + q1)(x^2 + rx + q2)(x^2 + sx + q3) ??? hard.
\\
\\ ALTERNATIVE: Use the Recillas correspondence.
\\ For a trigonal curve y^3 = f(x) of degree d, there's a correspondence
\\ to genus-2 curves. For j=0 elliptic curves, the relevant correspondence
\\ is through the theta-lift.

\\ 11. Try: the correct c is DETERMINED by the condition that
\\     y^2 = (x^3+b)(x^3+c) has char poly = (T^2-t*T+p)(T^2+t*T+p).
\\     This means: the TRACE of Jac(C) = 0 (biquadratic Weil poly).
\\ Condition: a1(Jac) = 0 AND a2(Jac) = 2p - t^2.
\\ For y^2 = (x^3+A)(x^3+B), the Frobenius trace a2 of Jac is related to
\\ the traces tA, tB of the constituent elliptic curves...
\\ Actually for Jac~E×E^t, a2 = 2p - tE^2.
\\
\\ The key formula from the theory of CM abelian surfaces:
\\ y^2 = (x^3+A)(x^3+B) has Jac ~ E1 x E2 (both j=0) iff
\\     A * B = (c1^3 + c2^3) for some specific (c1,c2)... unclear.

print("---- Algebraic check: what is 140*919 and 140*453 mod p ----");
print("140 * 453 mod p = ", lift(Mod(140*453, p)));
print("140 * 919 mod p = ", lift(Mod(140*919, p)));
print("453 * 919 mod p = ", lift(Mod(453*919, p)));
print("140 * 453 * 919 mod p = ", lift(Mod(140*453*919, p)));
print("");
print("b * b2 mod p = ", lift(Mod(b*b2, p)));
print("b * b2 / (140*453*919) mod p = ", lift(Mod(b*b2,p) * Mod(140*453*919,p)^(-1)));
print("");

\\ What is special about c=140?
print("c=140 factored: ", factor(140));
print("b=11 factored: ", factor(11));
print("b2=515 factored: ", factor(515));
print("d=11 factored: ", factor(11));
print("");
print("140 mod 3 = ", 140 % 3, "  b mod 3 = ", b % 3, "  b2 mod 3 = ", b2 % 3);
print("");

\\ Check if 140 = b * (something)^2 / d or similar
for(k=1, 100,
  for(m=1, p-1,
    if(lift(Mod(b * m^2, p)) == 140,
      print("140 = b * ", m, "^2 mod p = b * ", lift(Mod(m^2,p)));
      break(2)
    )
  )
);

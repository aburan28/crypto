\\ thread9_igusa_naive.gp
\\ Thread 9: Igusa invariants of naive-cover Jacobians y^2 = (x^3+b_i)(x^3+b_j)
\\
\\ For the "NONSPLIT-Q" simple Jacobians identified in thread8/9, compute
\\ Igusa-Clebsch invariants (I2:I4:I6:I10) using existing infrastructure.
\\
\\ Cases studied: those where a2-2p=-73 (CM field Q(sqrt(73), sqrt(-3)))
\\   p=19, pair (b1=2, b2=4): f(x) = x^6 + 6x^3 + 8
\\   p=37, pair (b1=2, b2=4): f(x) with g=2
\\   p=79, pair (b1=2, b2=4): f(x) with g=2
\\ Plus one control case: pair with different CM field.
\\
\\ The Igusa invariants (I2:I4:I6:I10) characterize Jac over Fbar_p.
\\ If two Jacobians at different p but same CM field give "CM-equivalent"
\\ invariants, that confirms the CM identification.
\\
\\ Run: gp -q thread9_igusa_naive.gp

default(parisize, 256000000);
default(timer, 0);

read("igusa_clebsch.gp");
read("igusa_clebsch_complete.gp");

\\ ========================================================
\\ Helper: Igusa invariants over F_p
\\ ========================================================

\\ Compute Clebsch A,B,C,D and Igusa I2,I4,I6,I10 for y^2=f(x) over F_p.
\\ f is given as a polynomial in x of degree 6 with F_p coefficients.
igusa_over_fp(f, p) = {
  my(A, B, C, D, I2, I4, I6, I10, res);
  A  = Mod(clebsch_A(f), p);
  B  = Mod(clebsch_B(f), p);
  C  = Mod(clebsch_C(f), p);
  D  = Mod(clebsch_D(f), p);
  I2  = Mod(-120, p) * A;
  I4  = Mod(-720, p) * A^2 + Mod(6750, p) * B;
  I6  = Mod(8640, p) * A^3 - Mod(108000, p) * A*B + Mod(202500, p) * C;
  I10 = Mod(-62208,p)*A^5 + Mod(972000,p)*A^3*B + Mod(1620000,p)*A^2*C
        - Mod(3037500,p)*A*B^2 - Mod(6075000,p)*B*C - Mod(4556250,p)*D;
  [lift(I2), lift(I4), lift(I6), lift(I10)]
};

\\ ========================================================
\\ Primitive root (smallest g with ord(g)=p-1)
\\ ========================================================
prim_root(p) = {
  my(phi, g, ok, facs, n);
  phi = p - 1;
  facs = factor(phi)[,1];
  for(g = 2, p-1,
    ok = 1;
    for(k = 1, #facs,
      if(Mod(g, p)^(phi/facs[k]) == 1, ok = 0; break)
    );
    if(ok, return(g))
  );
  error("no primitive root found")
};

\\ ========================================================
\\ Naive cover curve f for pair (i,j) of sextic twists
\\ ========================================================
naive_cover_f(p, i, j) = {
  my(g, bi, bj);
  g = prim_root(p);
  bi = Mod(g, p)^i;
  bj = Mod(g, p)^j;
  bi = lift(bi);
  bj = lift(bj);
  \\ f(x) = (x^3+bi)(x^3+bj) = x^6 + (bi+bj)*x^3 + bi*bj
  my(x);
  x = Mod('x, p);
  lift((x^3 + bi) * (x^3 + bj))
};

\\ ========================================================
\\ Main computation
\\ ========================================================

print("=============================================================");
print("Thread 9: Igusa invariants for naive-cover simple Jacobians");
print("=============================================================");
print();

\\ CM field Q(sqrt(73), sqrt(-3)) cases: a2-2p = -73
\\ Pairs (2,3) and (5,6) in 1-indexed = (1,2) and (4,5) in 0-indexed.
\\ With g as prim root: pair 0-indexed (1,2) means twist indices i=1,j=2.

cases_cm73 = [[19,1,2], [37,1,2], [79,1,2], [109,1,2]];

print("--- CM field K=Q(sqrt(73),sqrt(-3)) cases: a2-2p=-73 ---");
print("(Pairs with a1=0, D=4*73=292, sf=73)");
print();

for(k = 1, #cases_cm73,
  my(p, i, j, f, inv, I2, I4, I6, I10);
  p = cases_cm73[k][1];
  i = cases_cm73[k][2];
  j = cases_cm73[k][3];
  f = naive_cover_f(p, i, j);
  print("p=", p, " pair (", i, ",", j, ") [0-indexed]");
  print("  f(x) = ", f);
  inv = igusa_over_fp(f, p);
  I2 = inv[1]; I4 = inv[2]; I6 = inv[3]; I10 = inv[4];
  print("  Igusa: I2=", I2, "  I4=", I4, "  I6=", I6, "  I10=", I10);
  if(I10 != 0,
    \\ Normalize to projective ratios [I2:I4:I6:I10]
    \\ Absolute Igusa invariants: j1 = I2^5/I10, j2 = I2^3*I4/I10, j3 = I2^2*I6/I10
    my(j1, j2, j3);
    j1 = Mod(I2, p)^5 / Mod(I10, p);
    j2 = Mod(I2, p)^3 * Mod(I4, p) / Mod(I10, p);
    j3 = Mod(I2, p)^2 * Mod(I6, p) / Mod(I10, p);
    print("  j1=I2^5/I10 = ", lift(j1));
    print("  j2=I2^3*I4/I10 = ", lift(j2));
    print("  j3=I2^2*I6/I10 = ", lift(j3));
  ,
    print("  I10=0: curve has bad reduction or is not genus-2 smooth")
  );
  print();
);

\\ Non-CM-73 cases for comparison
print("--- Comparison cases (other CM fields) ---");
print();

\\ p=67, pair (1,2) [0-indexed]: a2-2p = 64-67*2... let me use the right pair
\\ From thread9 data: p=67, pair (1,2) has sf=265. Use pair (0,1) [0-indexed].
cases_other = [[31,0,3], [43,0,3], [61,0,3]];

for(k = 1, #cases_other,
  my(p, i, j, f, inv, I2, I4, I6, I10);
  p = cases_other[k][1];
  i = cases_other[k][2];
  j = cases_other[k][3];
  f = naive_cover_f(p, i, j);
  print("p=", p, " pair (", i, ",", j, ") [0-indexed]");
  print("  f(x) = ", f);
  inv = igusa_over_fp(f, p);
  I2 = inv[1]; I4 = inv[2]; I6 = inv[3]; I10 = inv[4];
  print("  Igusa: I2=", I2, "  I4=", I4, "  I6=", I6, "  I10=", I10);
  if(I10 != 0,
    my(j1, j2, j3);
    j1 = Mod(I2, p)^5 / Mod(I10, p);
    j2 = Mod(I2, p)^3 * Mod(I4, p) / Mod(I10, p);
    j3 = Mod(I2, p)^2 * Mod(I6, p) / Mod(I10, p);
    print("  j1=I2^5/I10 = ", lift(j1));
    print("  j2=I2^3*I4/I10 = ", lift(j2));
    print("  j3=I2^2*I6/I10 = ", lift(j3));
  ,
    print("  I10=0: curve has bad reduction or is not genus-2 smooth")
  );
  print();
);

\\ Check if I10 != 0 in all cm73 cases (I10=discriminant of curve, nonzero = smooth)
print("--- Smoothness check (I10 = discriminant proxy) ---");
for(k = 1, #cases_cm73,
  my(p, i, j, f, inv);
  p = cases_cm73[k][1];
  i = cases_cm73[k][2];
  j = cases_cm73[k][3];
  f = naive_cover_f(p, i, j);
  inv = igusa_over_fp(f, p);
  print("p=", p, " pair(", i, ",", j, "): I10=", inv[4], " (", if(inv[4]!=0, "smooth", "singular"), ")");
);
print();

print("Done.");

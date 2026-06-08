\\ howe_richelot_p43_v2.gp
\\ Rewrite of howe_richelot_p43.gp: avoid nested multi-line ff3_add() calls
\\ that fail in PARI 2.15.4.  Every intermediate value stored in a variable.
\\ Run: gp -q howe_richelot_p43_v2.gp

default(parisize, 64000000); default(timer, 0);

p = 43; z3 = 6; z3sq = 36;

\\ F_{43^3} = F_43[a]/(a^3 - 36).  Elements stored as [c0,c1,c2] = c0+c1*a+c2*a^2.

ff3_add(u,v) = [(u[1]+v[1])%p, (u[2]+v[2])%p, (u[3]+v[3])%p];
ff3_neg(u) = [(-u[1])%p, (-u[2])%p, (-u[3])%p];
ff3_scl(c,u) = [(c*u[1])%p, (c*u[2])%p, (c*u[3])%p];
ff3_c(c)    = [c%p, 0, 0];

ff3_mul(u,v) = {
  my(c0,c1,c2,c3,c4);
  c0 = (u[1]*v[1]) % p;
  c1 = (u[1]*v[2] + u[2]*v[1]) % p;
  c2 = (u[1]*v[3] + u[2]*v[2] + u[3]*v[1]) % p;
  c3 = (u[2]*v[3] + u[3]*v[2]) % p;
  c4 = (u[3]*v[3]) % p;
  [(c0 + 36*c3) % p,  (c1 + 36*c4) % p,  c2]
};

ff3_inv(u) = {
  my(a,b,c,N,ni);
  a=u[1]; b=u[2]; c=u[3];
  N = lift(Mod(a^3 + 36*b^3 + 36^2*c^3 - 3*36*a*b*c, p));
  if(N==0, error("ff3_inv: zero norm"));
  ni = lift(Mod(N,p)^(-1));
  [(a^2 - 36*b*c)*ni % p,  (36*c^2 - a*b)*ni % p,  (b^2 - a*c)*ni % p]
};

\\ Compute Det3 of the 3x3 matrix
\\   | 1    Gx1  Gc1 |
\\   | 1    Gx2  Gc2 |
\\   | 1    Gx3  Gc3 |
\\ (leading coeff of each G_i is 1 = [1,0,0])
\\ = (Gx2*Gc3 - Gx3*Gc2) - Gx1*(Gc3 - Gc2) + Gc1*(Gx3 - Gx2)
det3_G(Gx1,Gc1, Gx2,Gc2, Gx3,Gc3) = {
  my(t1,t2,t3,t4,t5,t6,t7,t8,t9,D);
  t1 = ff3_mul(Gx2, Gc3);
  t2 = ff3_mul(Gx3, Gc2);
  t3 = ff3_add(t1, ff3_neg(t2));
  t4 = ff3_add(Gc3, ff3_neg(Gc2));
  t5 = ff3_mul(Gx1, t4);
  t6 = ff3_add(t3, ff3_neg(t5));
  t7 = ff3_add(Gx3, ff3_neg(Gx2));
  t8 = ff3_mul(Gc1, t7);
  D = ff3_add(t6, t8);
  D
};

\\ H_num(j,k) = G_j' * G_k - G_j * G_k'  (numerator of Richelot H_i)
\\ G_i(x) = x^2 + Gxi*x + Gci  =>  G_i'(x) = 2x + Gxi
\\ G_j'*G_k: degree 3 (x^3 term = 2*1 = 2, but we only need x^2, x^1, x^0)
\\ G_j*G_k': degree 3 similarly.
\\ Difference (G_j'*G_k - G_j*G_k') has x^3 coeff = 0, so degree 2.
\\ x^2 coeff: 2*Gkx + Gjx - (Gjx + 2*Gkx) ... wait, more carefully:
\\   G_j'*G_k: x^2 coeff = 2*Gkx + Gjx*1
\\   G_j*G_k': x^2 coeff = 1*2 + Gjx*Gkx... no wait.
\\ Let me be explicit:
\\ G_j' = 2x + Gjx; G_k = x^2 + Gkx*x + Gkc
\\ G_j'*G_k = 2x^3 + (2*Gkx + Gjx)*x^2 + (2*Gkc + Gjx*Gkx)*x + Gjx*Gkc
\\ G_j = x^2 + Gjx*x + Gjc; G_k' = 2x + Gkx
\\ G_j*G_k' = 2x^3 + (Gkx + 2*Gjx)*x^2 + (Gjx*Gkx + 2*Gjc)*x + Gjc*Gkx
\\ Difference (H numerator):
\\ x^3: 0
\\ x^2: (2*Gkx + Gjx) - (Gkx + 2*Gjx) = Gkx - Gjx
\\ x^1: (2*Gkc + Gjx*Gkx) - (Gjx*Gkx + 2*Gjc) = 2*(Gkc - Gjc)
\\ x^0: Gjx*Gkc - Gjc*Gkx
H_num_coeffs(Gjx,Gjc, Gkx,Gkc) = {
  my(hx2, hx1, hx0);
  hx2 = ff3_add(Gkx, ff3_neg(Gjx));
  hx1 = ff3_scl(2, ff3_add(Gkc, ff3_neg(Gjc)));
  hx0 = ff3_add(ff3_mul(Gjx,Gkc), ff3_neg(ff3_mul(Gjc,Gkx)));
  [hx0, hx1, hx2]
};

\\ Multiply two degree-2 polys (each as [a0,a1,a2]) -> degree-4 poly [b0..b4]
poly2_mul(a, b) = {
  my(c0,c1,c2,c3,c4);
  c0 = ff3_mul(a[1],b[1]);
  c1 = ff3_add(ff3_mul(a[1],b[2]), ff3_mul(a[2],b[1]));
  c2 = ff3_add(ff3_add(ff3_mul(a[1],b[3]),ff3_mul(a[2],b[2])),ff3_mul(a[3],b[1]));
  c3 = ff3_add(ff3_mul(a[2],b[3]), ff3_mul(a[3],b[2]));
  c4 = ff3_mul(a[3],b[3]);
  [c0,c1,c2,c3,c4]
};

\\ Multiply degree-4 poly [a0..a4] by degree-2 poly [b0,b1,b2] -> degree-6 [c0..c6]
poly4_2_mul(a, b) = {
  my(c0,c1,c2,c3,c4,c5,c6);
  c0 = ff3_mul(a[1],b[1]);
  c1 = ff3_add(ff3_mul(a[1],b[2]), ff3_mul(a[2],b[1]));
  c2 = ff3_add(ff3_add(ff3_mul(a[1],b[3]),ff3_mul(a[2],b[2])),ff3_mul(a[3],b[1]));
  c3 = ff3_add(ff3_add(ff3_mul(a[2],b[3]),ff3_mul(a[3],b[2])),ff3_mul(a[4],b[1]));
  c4 = ff3_add(ff3_add(ff3_mul(a[3],b[3]),ff3_mul(a[4],b[2])),ff3_mul(a[5],b[1]));
  c5 = ff3_add(ff3_mul(a[4],b[3]),ff3_mul(a[5],b[2]));
  c6 = ff3_mul(a[5],b[3]);
  [c0,c1,c2,c3,c4,c5,c6]
};

\\ Given three quadratics [Gix_c, Gix_x, Gix_x2] with Gix_x2=[1,0,0] (monic),
\\ compute H1*H2*H3 over F_{43^3} and check it's in F_43[x].
\\ Returns [ok, a3, b] where y^2 = x^6 + a3*x^3 + b (monic Z/3Z form),
\\ or [-1, ...] on failure.
richelot_sigma(Gx1,Gc1, Gx2,Gc2, Gx3,Gc3) = {
  my(Delta, Di);
  my(H1n, H2n, H3n, H1, H2, H3);
  my(H1H2, H1H2H3);
  my(ok, lc, lci, a3, b);
  Delta = det3_G(Gx1,Gc1, Gx2,Gc2, Gx3,Gc3);
  if(Delta[2]!=0 || Delta[3]!=0,
    print("Delta not in F_43: ", Delta);
    return([-1,0,0])
  );
  Di = ff3_inv(Delta);
  H1n = H_num_coeffs(Gx2,Gc2, Gx3,Gc3);   \\ cyclic: H1 uses (j=2,k=3)
  H2n = H_num_coeffs(Gx3,Gc3, Gx1,Gc1);   \\ H2 uses (j=3,k=1)
  H3n = H_num_coeffs(Gx1,Gc1, Gx2,Gc2);   \\ H3 uses (j=1,k=2)
  H1 = [ff3_mul(H1n[1],Di), ff3_mul(H1n[2],Di), ff3_mul(H1n[3],Di)];
  H2 = [ff3_mul(H2n[1],Di), ff3_mul(H2n[2],Di), ff3_mul(H2n[3],Di)];
  H3 = [ff3_mul(H3n[1],Di), ff3_mul(H3n[2],Di), ff3_mul(H3n[3],Di)];
  H1H2 = poly2_mul(H1, H2);
  H1H2H3 = poly4_2_mul(H1H2, H3);
  \\ Check all coefficients are in F_43
  ok = 1;
  for(i=1, 7,
    if(H1H2H3[i][2]!=0 || H1H2H3[i][3]!=0,
      print("  coeff x^",i-1," not in F_43: ", H1H2H3[i]);
      ok = 0
    )
  );
  if(!ok, return([-1,0,0]));
  lc = H1H2H3[7][1];
  if(lc==0, print("  leading coeff=0"); return([-1,0,0]));
  lci = lift(Mod(lc,p)^(-1));
  a3 = (H1H2H3[4][1]*lci) % p;
  b  = (H1H2H3[1][1]*lci) % p;
  [1, a3, b]
};

\\ ================================================================
print("================================================================");
print("Howe-glued genus-2 curve over F_43: Richelot (2,2)-isogeny");
print("E1: y^2=x^3+7  (trace 13)   E2: y^2=x^3+13 (trace -13)");
print("beta = 2*alpha  (since (2a)^3=8*(-7)=56=30=-13 mod 43)");
print("================================================================");
print("");

\\ Verify basic setup
E1 = ellinit([0,7],p);  E2 = ellinit([0,13],p);
t1 = p+1-ellcard(E1);   t2 = p+1-ellcard(E2);
print("E1 trace=",t1,"  E2 trace=",t2,"  (expected +/-13)");
print("(2*alpha)^3 mod 43 = ",(8*36)%p,"  (expected 30=-13)");
print("");

\\ Basis: alpha = [0,1,0] in F_{43^3}=F_43[a]/(a^3-36)
alpha  = [0,1,0];
alpha2 = [0,0,1];

\\ ================================================================
\\ sigma_0: pair (alpha, 2*alpha), (z3*alpha, z3*2*alpha), (z3^2*alpha, z3^2*2*alpha)
\\ G_i = (x - z3^{i-1}*alpha)(x - z3^{i-1}*2*alpha)
\\     = x^2 - (3*z3^{i-1}*alpha)*x + (2*z3^{2(i-1)}*alpha^2)
\\ sigma_0: s0=3*alpha, q0=2*alpha^2
\\ ================================================================

print("---- sigma_0 ----");
s0 = ff3_scl(3, alpha);     \\ 3*alpha
q0 = ff3_scl(2, alpha2);    \\ 2*alpha^2

G0_1x = ff3_neg(s0);               G0_1c = q0;
G0_2x = ff3_neg(ff3_scl(z3,s0));   G0_2c = ff3_scl(z3sq,q0);
G0_3x = ff3_neg(ff3_scl(z3sq,s0)); G0_3c = ff3_scl(z3,q0);

D0 = det3_G(G0_1x,G0_1c, G0_2x,G0_2c, G0_3x,G0_3c);
print("Delta(sigma_0) = ", D0, "  (expected [39,0,0])");

r0 = richelot_sigma(G0_1x,G0_1c, G0_2x,G0_2c, G0_3x,G0_3c);
if(r0[1]==1,
  print("sigma_0 Richelot image: y^2 = x^6 + ", r0[2], "*x^3 + ", r0[3]);
  cp0 = hyperellcharpoly((x^6 + r0[2]*x^3 + r0[3]) * Mod(1,p));
  print("  char poly = ", cp0);
  print("  target    = T^4-83*T^2+1849")
,
  print("sigma_0 FAILED")
);
print("");

\\ ================================================================
\\ sigma_1: pair (alpha, z3*2*alpha), (z3*alpha, z3^2*2*alpha), (z3^2*alpha, 2*alpha)
\\ s1 = alpha + z3*2*alpha = (1+2*z3)*alpha = 13*alpha
\\ q1 = alpha * z3*2*alpha = 2*z3*alpha^2 = 12*alpha^2
\\ ================================================================

print("---- sigma_1 ----");
s1 = ff3_scl((1+2*z3)%p, alpha);    \\ 13*alpha
q1 = ff3_scl((2*z3)%p, alpha2);     \\ 12*alpha^2

G1_1x = ff3_neg(s1);               G1_1c = q1;
G1_2x = ff3_neg(ff3_scl(z3,s1));   G1_2c = ff3_scl(z3sq,q1);
G1_3x = ff3_neg(ff3_scl(z3sq,s1)); G1_3c = ff3_scl(z3,q1);

D1 = det3_G(G1_1x,G1_1c, G1_2x,G1_2c, G1_3x,G1_3c);
print("Delta(sigma_1) = ", D1, "  (expected [37,0,0])");

r1 = richelot_sigma(G1_1x,G1_1c, G1_2x,G1_2c, G1_3x,G1_3c);
if(r1[1]==1,
  print("sigma_1 Richelot image: y^2 = x^6 + ", r1[2], "*x^3 + ", r1[3]);
  cp1 = hyperellcharpoly((x^6 + r1[2]*x^3 + r1[3]) * Mod(1,p));
  print("  char poly = ", cp1);
  print("  target    = T^4-83*T^2+1849")
,
  print("sigma_1 FAILED")
);
print("");

\\ ================================================================
\\ sigma_2: pair (alpha, z3^2*2*alpha), (z3*alpha, 2*alpha), (z3^2*alpha, z3*2*alpha)
\\ s2 = alpha + z3^2*2*alpha = (1+2*z3^2)*alpha = (1+72)*alpha = 73*alpha = 30*alpha
\\ q2 = alpha * z3^2*2*alpha = 2*z3^2*alpha^2 = 72*alpha^2 = 29*alpha^2
\\ ================================================================

print("---- sigma_2 ----");
s2 = ff3_scl((1+2*z3sq)%p, alpha);    \\ (1+72)*alpha = 30*alpha
q2 = ff3_scl((2*z3sq)%p, alpha2);     \\ 72*alpha^2 = 29*alpha^2

G2_1x = ff3_neg(s2);               G2_1c = q2;
G2_2x = ff3_neg(ff3_scl(z3,s2));   G2_2c = ff3_scl(z3sq,q2);
G2_3x = ff3_neg(ff3_scl(z3sq,s2)); G2_3c = ff3_scl(z3,q2);

D2 = det3_G(G2_1x,G2_1c, G2_2x,G2_2c, G2_3x,G2_3c);
print("Delta(sigma_2) = ", D2);

r2 = richelot_sigma(G2_1x,G2_1c, G2_2x,G2_2c, G2_3x,G2_3c);
if(r2[1]==1,
  print("sigma_2 Richelot image: y^2 = x^6 + ", r2[2], "*x^3 + ", r2[3]);
  cp2 = hyperellcharpoly((x^6 + r2[2]*x^3 + r2[3]) * Mod(1,p));
  print("  char poly = ", cp2);
  print("  target    = T^4-83*T^2+1849")
,
  print("sigma_2 FAILED")
);
print("");

\\ ================================================================
\\ Match results against the 7 canonical Z/3Z classes
\\ ================================================================

print("================================================================");
print("7 canonical Z/3Z classes (from 2026-06-01 search):");
print("================================================================");
classes_a = [1, 2, 3, 4, 6, 8, 16];
classes_b = [2, 8, 42, 32, 39, 42, 39];

check_match(label, res) = {
  my(found, ca, cb);
  if(res[1]!=1, print(label,": computation failed"); return());
  found = 0;
  for(i=1,7,
    ca=classes_a[i]; cb=classes_b[i];
    if(res[2]==ca && res[3]==cb,
      print(label, ": MATCH class ", i, " (a=", ca, ", b=", cb, ") ***");
      found = 1
    )
  );
  if(!found,
    print(label, ": a=", res[2], " b=", res[3], " — NOT in canonical 7 classes")
  )
};

check_match("sigma_0", r0);
check_match("sigma_1", r1);
check_match("sigma_2", r2);
print("");

\\ ================================================================
\\ Verify the matched curve's Igusa invariants against the table
\\ ================================================================
print("================================================================");
print("Igusa-Clebsch invariants (from igusa_7classes_full.gp table):");
print("  Class | a  | b  | I2 | I4 | I6 | I10");
print("  1     | 1  | 2  | 21 | 33 | 41 | 35");
print("  2     | 2  | 8  | 41 | 12 | 1  | 21");
print("  3     | 3  | 42 | 18 | 37 | 22 | 35");
print("  4     | 4  | 32 | 35 | 20 | 21 | 4 ");
print("  5     | 6  | 39 | 29 | 33 | 32 | 21");
print("  6     | 8  | 42 | 11 | 19 | 11 | 11");
print("  7     | 16 | 39 | 1  | 3  | 16 | 41");
print("================================================================");
print("");
print("Summary:");
if(r0[1]==1, print("  sigma_0 -> (a=",r0[2],", b=",r0[3],")"));
if(r1[1]==1, print("  sigma_1 -> (a=",r1[2],", b=",r1[3],")"));
if(r2[1]==1, print("  sigma_2 -> (a=",r2[2],", b=",r2[3],")"));
print("");
print("================================================================");
print("Done.");
print("================================================================");

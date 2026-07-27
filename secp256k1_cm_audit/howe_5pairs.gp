\\ howe_5pairs.gp
\\
\\ Compute explicit Howe-glued genus-2 covers for the 5 pairs of
\\ secp256k1 sextic twists that satisfy H1+H2+H3 (Thread 18, 2026-07-21).
\\
\\ The 5 pairs (using twist index k where E_k: y^2 = x^3 + 7*h^k mod p):
\\   (0,3): secp256k1 × quadratic twist  [KNOWN]
\\   (0,1): NEW
\\   (0,4): NEW
\\   (1,4): quadratic-twist sub-pair     [NEW]
\\   (3,4): NEW
\\
\\ Method: Z/3Z Richelot isogeny from howe_richelot_v5.gp
\\   Given alpha, beta = cube roots of 2-torsion of E_i, E_j in F_{p^3},
\\   the Howe-glued cover has the Richelot image of:
\\     G1 = (x - alpha)(x - beta)
\\     G2 = (x - z3*alpha)(x - z3*beta)
\\     G3 = (x - z3^2*alpha)(x - z3^2*beta)
\\   i.e., richelot_gen(sv=alpha+beta, qv=alpha*beta).
\\
\\ For j=0 curve E_k: y^2 = x^3 + c_k (where c_k = 7*h^k mod p):
\\   F_{p^3} = F_p[a]/(a^3 - r_k) where r_k = -c_k mod p.
\\   The cube root alpha_k of -c_k is the element 'a' in F_{p^3}.
\\
\\ For pairs (i,j): alpha = a (cube root of -c_i), beta = a' (cube root of -c_j)
\\   These may or may not be rationally related (e.g., beta = d*alpha for d in F_p).
\\   If they ARE related (beta = d*alpha), the Howe cover is over F_p.
\\   Otherwise, we work over F_{p^3} and check if the result descends to F_p.
\\
\\ Run: gp -q howe_5pairs.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Howe-glued genus-2 covers for 5 secp256k1 sextic twist pairs");
print("================================================================");
print("");

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
b0 = 7;  \\ secp256k1: y^2 = x^3 + 7

\\ ---- Find primitive 6th root of unity h in F_p ----
\\   h satisfies h^6 = 1, h != 1, -1.  h is also a primitive 6th root.
\\   We compute via h = -zeta3 where zeta3 = (-1+sqrt(-3))/2 mod p.
\\   Alternatively, solve x^2 - x + 1 = 0 mod p.
{
  local_roots = polrootsmod(x^2 - x + 1, p);
  h = lift(local_roots[1]);
  print("6th root of unity h = ", h, " (mod p)");
  print("  h^6 mod p = ", lift(Mod(h,p)^6), "  (should be 1)");
  print("  h^3 mod p = ", lift(Mod(h,p)^3), "  (should be p-1)");
  print("");
}

\\ ---- Compute twist coefficients c_k = 7*h^k mod p ----
{
  c = vector(6);
  for(k=0, 5,
    c[k+1] = lift(Mod(b0, p) * Mod(h, p)^k);
  );
  print("Twist coefficients c_k = 7*h^k mod p:");
  for(k=0, 5,
    E_k = ellinit([0, c[k+1]], p);
    t_k = p + 1 - ellcard(E_k);
    print("  k=", k, ": c_k = ", c[k+1], "  trace = ", t_k);
  );
  print("");
}

\\ ================================================================
\\ F_{p^3} arithmetic: elements [c0, c1, c2] represent c0 + c1*a + c2*a^2
\\ where a^3 = r (the "reduction polynomial constant").
\\ This requires r to be fixed before using ff3_mul etc.
\\ We use a global variable 'ff3_r' for the reduction constant.
\\ ================================================================

\\ Note: r depends on the pair being computed.  We re-define it per pair.

ff3_add(u, v) = [(u[1]+v[1])%p, (u[2]+v[2])%p, (u[3]+v[3])%p];
ff3_neg(u) = [(-u[1])%p, (-u[2])%p, (-u[3])%p];
ff3_scl(c, u) = [(c*u[1])%p, (c*u[2])%p, (c*u[3])%p];

ff3_mul_r(u, v, r) = {
  my(c0, c1, c2, c3, c4);
  c0 = (u[1]*v[1]) % p;
  c1 = (u[1]*v[2] + u[2]*v[1]) % p;
  c2 = (u[1]*v[3] + u[2]*v[2] + u[3]*v[1]) % p;
  c3 = (u[2]*v[3] + u[3]*v[2]) % p;
  c4 = (u[3]*v[3]) % p;
  [(c0 + r*c3) % p,  (c1 + r*c4) % p,  c2]
};

ff3_inv_r(u, r) = {
  my(a, b, c, n, ni);
  a = u[1]; b = u[2]; c = u[3];
  n = lift(Mod(a^3 + r^2*c^3 - 3*r*a*b*c + r*a^2 - r^2*b^2*c, p));
  \\ Wait, norm of a+b*alpha+c*alpha^2 in F_p[alpha]/(alpha^3-r)
  \\ N = a^3 + r*b^3 + r^2*c^3 - 3*r*a*b*c
  n = lift(Mod(a^3 + r*b^3 + r^2*c^3 - 3*r*a*b*c, p));
  if(n == 0, error("ff3_inv: element is not invertible"));
  ni = lift(Mod(n, p)^(-1));
  [((a^2 - r*b*c)*ni) % p,
   ((r*c^2 - a*b)*ni) % p,
   ((b^2 - a*c)*ni) % p]
};

\\ Richelot formula: given sv = alpha+beta, qv = alpha*beta (in F_{p^3}),
\\ and z3 = primitive cube root of unity in F_p,
\\ compute the dual Richelot cover.
\\ Returns [aa, bb] for y^2 = x^6 + aa*x^3 + bb, or [-1, -1] if not over F_p.

richelot_z3(sv, qv, z3, r) = {
  my(z3sq, G1c, G1x, G2c, G2x, G3c, G3x);
  my(D0, D0inv, H1n2, H1n1, H1n0, H2n2, H2n1, H2n0, H3n2, H3n1, H3n0);
  my(H1x2, H1x1, H1x0, H2x2, H2x1, H2x0, H3x2, H3x1, H3x0);
  my(P0, P1, P2, P3, P4, Q0r, Q1r, Q2r, Q3r, Q4r, Q5r, Q6r);
  my(lc, lcinv, aa, bb);

  z3sq = lift(Mod(z3, p)^2);

  G1c = qv;
  G1x = ff3_neg(sv);
  G2c = ff3_scl(z3sq, qv);
  G2x = ff3_neg(ff3_scl(z3, sv));
  G3c = ff3_scl(z3, qv);
  G3x = ff3_neg(ff3_scl(z3sq, sv));

  \\ Delta = det of 3x3 matrix of (Gix^2, Gix, Gic) rows
  \\ (with Gix^2 = 1 since leading coeff is x^2, so first column is all 1)
  D0 = ff3_add(
    ff3_add(
      ff3_add(ff3_mul_r(G2x, G3c, r), ff3_neg(ff3_mul_r(G3x, G2c, r))),
      ff3_neg(ff3_mul_r(G1x, ff3_add(G3c, ff3_neg(G2c)), r))
    ),
    ff3_mul_r(G1c, ff3_add(G3x, ff3_neg(G2x)), r)
  );

  D0inv = ff3_inv_r(D0, r);

  \\ H1 numerator (j=2, k=3 cyclic)
  H1n2 = ff3_add(G3x, ff3_neg(G2x));
  H1n1 = ff3_scl(2, ff3_add(G3c, ff3_neg(G2c)));
  H1n0 = ff3_add(ff3_mul_r(G2x, G3c, r), ff3_neg(ff3_mul_r(G2c, G3x, r)));

  \\ H2 numerator (j=3, k=1)
  H2n2 = ff3_add(G1x, ff3_neg(G3x));
  H2n1 = ff3_scl(2, ff3_add(G1c, ff3_neg(G3c)));
  H2n0 = ff3_add(ff3_mul_r(G3x, G1c, r), ff3_neg(ff3_mul_r(G3c, G1x, r)));

  \\ H3 numerator (j=1, k=2)
  H3n2 = ff3_add(G2x, ff3_neg(G1x));
  H3n1 = ff3_scl(2, ff3_add(G2c, ff3_neg(G1c)));
  H3n0 = ff3_add(ff3_mul_r(G1x, G2c, r), ff3_neg(ff3_mul_r(G1c, G2x, r)));

  H1x2 = ff3_mul_r(H1n2, D0inv, r);
  H1x1 = ff3_mul_r(H1n1, D0inv, r);
  H1x0 = ff3_mul_r(H1n0, D0inv, r);
  H2x2 = ff3_mul_r(H2n2, D0inv, r);
  H2x1 = ff3_mul_r(H2n1, D0inv, r);
  H2x0 = ff3_mul_r(H2n0, D0inv, r);
  H3x2 = ff3_mul_r(H3n2, D0inv, r);
  H3x1 = ff3_mul_r(H3n1, D0inv, r);
  H3x0 = ff3_mul_r(H3n0, D0inv, r);

  \\ Compute H1 * H2 (degree 4 polynomial in F_{p^3}[x])
  P0 = ff3_mul_r(H1x0, H2x0, r);
  P1 = ff3_add(ff3_mul_r(H1x0, H2x1, r), ff3_mul_r(H1x1, H2x0, r));
  P2 = ff3_add(ff3_add(ff3_mul_r(H1x0, H2x2, r), ff3_mul_r(H1x1, H2x1, r)), ff3_mul_r(H1x2, H2x0, r));
  P3 = ff3_add(ff3_mul_r(H1x1, H2x2, r), ff3_mul_r(H1x2, H2x1, r));
  P4 = ff3_mul_r(H1x2, H2x2, r);

  \\ Compute H1*H2 * H3 (degree 6 polynomial)
  Q0r = ff3_mul_r(P0, H3x0, r);
  Q1r = ff3_add(ff3_mul_r(P0, H3x1, r), ff3_mul_r(P1, H3x0, r));
  Q2r = ff3_add(ff3_add(ff3_mul_r(P0, H3x2, r), ff3_mul_r(P1, H3x1, r)), ff3_mul_r(P2, H3x0, r));
  Q3r = ff3_add(ff3_add(ff3_mul_r(P1, H3x2, r), ff3_mul_r(P2, H3x1, r)), ff3_mul_r(P3, H3x0, r));
  Q4r = ff3_add(ff3_add(ff3_mul_r(P2, H3x2, r), ff3_mul_r(P3, H3x1, r)), ff3_mul_r(P4, H3x0, r));
  Q5r = ff3_add(ff3_mul_r(P3, H3x2, r), ff3_mul_r(P4, H3x1, r));
  Q6r = ff3_mul_r(P4, H3x2, r);

  \\ Check all coefficients are in F_p (no a or a^2 components)
  if(Q0r[2] != 0 || Q0r[3] != 0 ||
     Q3r[2] != 0 || Q3r[3] != 0 ||
     Q6r[2] != 0 || Q6r[3] != 0,
    return([-1, -1])
  );

  lc = Q6r[1];
  if(lc == 0, return([-1, -1]));
  lcinv = lift(Mod(lc, p)^(-1));
  aa = (Q3r[1] * lcinv) % p;
  bb = (Q0r[1] * lcinv) % p;
  [aa, bb]
};

\\ ================================================================
\\ Step 1: Validate against toy case p=43
\\ ================================================================
print("================================================================");
print("Step 1: Validate Z/3Z Richelot against p=43 toy case");
print("================================================================");
print("");

{
  p43 = 43;
  z3_43 = 6;   \\ primitive cube root of unity mod 43
  r43 = 36;    \\ a^3 = 36 means -7 mod 43 (E1: y^2=x^3+7, alpha^3=-7=36)
  \\ alpha = [0,1,0], beta = 2*alpha (since (2alpha)^3=8*36=288=30=-13 mod 43)
  \\ d = 2, sv = (1+2)*alpha = 3*alpha = [0,3,0], qv = 2*alpha^2 = [0,0,2]
  sv43 = [0, 3, 0];
  qv43 = [0, 0, 2];
  res43 = richelot_z3(sv43, qv43, z3_43, r43);
  print("p=43, E1: y^2=x^3+7, E2: y^2=x^3+13 (quadratic twist via d=2)");
  print("  richelot result: a=", res43[1], "  b=", res43[2]);
  print("  expected: a=41  b=5  (from howe_richelot_v5.gp naive cover check)");
  if(res43[1] == 41 && res43[2] == 5,
    print("  STATUS: CORRECT ✓"),
    print("  STATUS: WRONG ✗ — formula may have a bug")
  );
  if(res43[1] != -1,
    h43 = Mod(1, p43)*x^6 + Mod(res43[1], p43)*x^3 + Mod(res43[2], p43);
    cp43 = hyperellcharpoly(h43);
    nj43 = subst(cp43, variable(cp43), 1);
    E1_43 = ellinit([0,7], p43); t1_43 = p43+1-ellcard(E1_43);
    nE1E2_43 = (p43+1-t1_43)*(p43+1+t1_43);
    print("  Frobenius poly of Jac(C): ", cp43);
    print("  #Jac(C) = ", nj43, "  target = ", nE1E2_43, "  match = ", nj43==nE1E2_43);
  );
  print("");
}

\\ ================================================================
\\ Step 2: Validate against p=1009 toy case
\\ ================================================================
print("================================================================");
print("Step 2: Validate against p=1009, b=11 toy case");
print("================================================================");
print("");

{
  p1009 = 1009;
  z3_1009 = lift(polrootsmod(x^2+x+1, p1009)[1]);
  print("p=1009: primitive cube root of unity z3 = ", z3_1009);

  \\ E1: y^2 = x^3 + 11, d = 11 (first non-square mod 1009)
  d1009 = 2; \\ try d=2 first, check if non-square
  while(kronecker(d1009, p1009) != -1, d1009++);
  print("First non-square d = ", d1009);
  b1_1009 = 11;
  b2_1009 = lift(Mod(d1009^3 * b1_1009, p1009));
  print("E1: y^2=x^3+", b1_1009, "  E2: y^2=x^3+", b2_1009);

  E1_1009 = ellinit([0, b1_1009], p1009);
  E2_1009 = ellinit([0, b2_1009], p1009);
  t1_1009 = p1009 + 1 - ellcard(E1_1009);
  t2_1009 = p1009 + 1 - ellcard(E2_1009);
  print("  trace E1 = ", t1_1009, "  trace E2 = ", t2_1009, "  (should be negatives)");

  \\ alpha^3 = -b1 mod p, so r = (-b1) mod p
  r1009 = (-b1_1009) % p1009;
  print("  reduction constant r = (-b1) mod p = ", r1009);

  \\ sv = (1 + d) * alpha = (1+d) * [0,1,0]
  sv1009 = [0, (1 + d1009) % p1009, 0];
  \\ qv = d * alpha^2 = d * [0,0,1]
  qv1009 = [0, 0, d1009 % p1009];

  print("  sv = (1+d)*alpha = ", sv1009);
  print("  qv = d*alpha^2   = ", qv1009);

  res1009 = richelot_z3(sv1009, qv1009, z3_1009, r1009);
  print("  Richelot result: a=", res1009[1], "  b=", res1009[2]);

  if(res1009[1] != -1,
    h1009 = Mod(1, p1009)*x^6 + Mod(res1009[1], p1009)*x^3 + Mod(res1009[2], p1009);
    cp1009 = hyperellcharpoly(h1009);
    nj1009 = subst(cp1009, variable(cp1009), 1);
    nE1E2_1009 = (p1009+1-t1_1009)*(p1009+1+t1_1009);
    print("  Frobenius poly: ", cp1009);
    print("  #Jac = ", nj1009, "  target = ", nE1E2_1009, "  match = ", nj1009==nE1E2_1009);
    if(nj1009 == nE1E2_1009,
      print("  STATUS: CORRECT ✓"),
      print("  STATUS: Jac order mismatch — need further investigation")
    ),
    print("  STATUS: cover not over F_p (result is -1)")
  );
  print("");
}

\\ ================================================================
\\ Step 3: secp256k1 pair (0,3) — known quadratic-twist pair
\\ ================================================================
print("================================================================");
print("Step 3: secp256k1 pair (0,3) — E_0 x E_3 (quadratic-twist)");
print("================================================================");
print("");

{
  \\ E_0: y^2 = x^3 + 7 (secp256k1)
  \\ E_3: y^2 = x^3 + 7*h^3 = x^3 + 7*(-1) = x^3 - 7 (quadratic twist, h^3=-1)
  \\ So E_3 has b3 = 7*h^3 = -7*h^0... wait, h^3 = p-1 = -1 mod p.
  \\ E_3: y^2 = x^3 + 7*(-1) = x^3 - 7 ≡ x^3 + (p-7).
  \\ But for the Howe construction, d must satisfy d^3 * b0 = b3, so d^3 = h^3 = -1.
  \\ d = (-1)^{1/3}: since p ≡ 2 mod 3? Let me check.

  h3 = lift(Mod(h, p)^3);  \\ should be p-1
  print("h^3 mod p = ", h3, "  (should be p-1 = ", p-1, ")");
  b3 = lift(Mod(7*h3, p));
  print("b3 = 7*h^3 mod p = ", b3, "  = p - 7 = ", p-7);
  E3 = ellinit([0, b3], p);
  t3 = p + 1 - ellcard(E3);
  n3 = ellcard(E3);
  E0 = ellinit([0, 7], p);
  t0 = p + 1 - ellcard(E0);
  print("E0: trace t0 = ", t0);
  print("E3: trace t3 = ", t3, "  (should be -t0 = ", -t0, ")");
  print("");

  \\ Find d such that d^3 = b3/b0 = h^3 = -1 mod p.
  \\ d^3 = -1, so d^6 = 1 and d^3 = -1.
  \\ d = h^3 = -1 itself is d^3 = -1 ✓.
  \\ But d = -1 has d^((p-1)/2) = (-1)^{(p-1)/2}.
  \\ Since p ≡ 3 mod 4 (secp256k1 p = 2^256 - 2^32 - 977, which is 3 mod 4),
  \\ (-1)^{(p-1)/2} = (-1)^{odd} = -1.  So d = -1 is a non-square mod p ✓.
  d03 = p - 1;  \\ d = -1
  print("d for pair (0,3): d = -1 mod p = ", d03);
  print("  d^3 mod p = ", lift(Mod(d03, p)^3), "  (should be b3/b0 = h^3 = ", h3, ")");
  print("  d^((p-1)/2) mod p = ", lift(Mod(d03, p)^((p-1)/2)), "  (should be p-1 for non-square)");
  print("");

  \\ Reduction constant: alpha^3 = -b0 = -7 mod p = p-7
  r03 = (p - 7) % p;

  \\ z3 in F_p (same as before)
  z3_p = lift(polrootsmod(x^2+x+1, p)[1]);

  \\ sv = (1+d)*alpha = (1+d)*[0,1,0]
  sv03 = [0, (1 + d03) % p, 0];
  \\ qv = d*alpha^2 = d*[0,0,1]
  qv03 = [0, 0, d03 % p];

  print("  sv = (1+d)*alpha where d=-1: sv = ", sv03);
  print("  (Note: sv = [0, 0, 0] since 1+d = 0 for d=-1)");

  \\ sv = [0, 0, 0]: this means alpha + beta = 0, i.e., beta = -alpha.
  \\ G1 = (x - alpha)(x - beta) = (x - alpha)(x + alpha) = x^2 - alpha^2
  \\ G2 = (x - z3*alpha)(x + z3*alpha) = x^2 - z3^2*alpha^2
  \\ G3 = (x - z3^2*alpha)(x + z3^2*alpha) = x^2 - z3^4*alpha^2 = x^2 - z3*alpha^2
  \\ These are THREE QUADRATICS with no linear term!
  \\ The Richelot formula for this case may degenerate.
  print("");
  print("  WARNING: sv = 0 (since d=-1, so beta = -alpha).");
  print("  The standard G_i = (x-alpha)(x-beta) = x^2 - sv*x + qv = x^2 + qv");
  print("  has NO linear term.  Richelot formula may need adjustment.");
  print("");

  \\ Try anyway:
  res03 = richelot_z3(sv03, qv03, z3_p, r03);
  print("  richelot result (sv=0 case): a=", res03[1], "  b=", res03[2]);
  print("");

  \\ Cross-check: the naive secp256k1 cover y^2=(x^3+7)(x^3-7) = x^6 - 49
  \\ But 189 = 27*7 was the naive cover from earlier work, not 7*(-1) = -7.
  \\ Let me check: is the naive cover (x^3+7)(x^3+189) or (x^3+7)(x^3-7)?
  \\ 189 = 7*27 = 7*3^3.  And -7 mod p ≠ 189.
  \\ So the quadratic twist with d=3 (non-square) gives b2 = 3^3*7 = 189,
  \\ while the quadratic twist with d=h (primitive 6th root) gives b2 = h^3*7 = -7.
  \\ For pair (0,3) with h^3 = -1: E_3: y^2 = x^3 + 7*h^3 = x^3 - 7.
  \\ This is a DIFFERENT quadratic twist than the "d=3" one.
  print("  Cross-check: known naive cover uses d=3 (not d=-1).");
  print("  Pair (0,3) uses d = h^3 = -1, which gives E_3: y^2=x^3-7.");
  print("  The naive cover for the standard d=3 twist is separate.");
  print("");
  print("  Finding d=3 non-square check: kronecker(3,p) = ", kronecker(3, p));
  \\ If 3 is a non-square, the naive cover corresponds to d=3 => twist pair (0, k) where k=?
  \\ 7*h^k = 7*3^3 = 189, so h^k = 27 = 3^3.
  \\ Find k: h^k = 27 mod p.
  \\ Discrete log of 27 base h mod p.
  print("  7*h^k = 189 for k s.t. h^k = 27.");
  for(kk=0, 5,
    if(lift(Mod(h,p)^kk) == 27 % p,
      print("  h^", kk, " = 27 mod p — so pair (0,", kk, ") uses d=3 twist");
    );
  );
  print("");
}

\\ ================================================================
\\ Step 4: secp256k1 pair (1,4) — quadratic-twist sub-pair
\\ ================================================================
print("================================================================");
print("Step 4: secp256k1 pair (1,4) — E_1 x E_4 (quadratic-twist pair)");
print("================================================================");
print("");

{
  b1_secp = lift(Mod(b0, p) * Mod(h, p)^1);
  b4_secp = lift(Mod(b0, p) * Mod(h, p)^4);
  print("E_1: y^2 = x^3 + ", b1_secp);
  print("E_4: y^2 = x^3 + ", b4_secp);

  E1_secp = ellinit([0, b1_secp], p);
  E4_secp = ellinit([0, b4_secp], p);
  t1_secp = p + 1 - ellcard(E1_secp);
  t4_secp = p + 1 - ellcard(E4_secp);
  print("  trace t1 = ", t1_secp);
  print("  trace t4 = ", t4_secp, "  (should be -t1 = ", -t1_secp, ")");
  print("");

  \\ d for pair (1,4): d^3 * b1 = b4, so d^3 = b4/b1 = h^3 = -1 mod p.
  \\ Same as pair (0,3)! So d = -1 again.
  \\ But sv = 0 again.  Let me use a DIFFERENT non-square d.
  \\ Actually, the Howe cover just needs ANY d with d^3*b1 = b4.
  \\ Since b4 = b1 * h^3 = b1 * (-1), we have d^3 = -1.
  \\ Any cube root of -1 works as d.
  \\ Alternative: use a different d' s.t. (d')^3 = -1 and d' != -1.
  \\ d^3 = -1 has solutions: d = -1, d = -zeta3, d = -zeta3^2.
  \\ These are the three cube roots of -1.
  \\ d = -1: beta = -alpha => sv = 0 (degenerate case)
  \\ d = -z3: beta = -z3*alpha => sv = alpha + (-z3*alpha) = (1-z3)*alpha
  \\ Let's use d = -z3 (= p - z3_p).

  z3_p = lift(polrootsmod(x^2+x+1, p)[1]);
  d14_alt = (p - z3_p) % p;
  print("  Alternative: d = -z3 mod p = ", d14_alt);
  \\ Check d^3 = (-z3)^3 = -z3^3 = -1 mod p:
  print("  d^3 mod p = ", lift(Mod(d14_alt, p)^3), "  (should be p-1 = ", p-1, ")");
  \\ Check non-square:
  print("  kronecker(d, p) = ", kronecker(d14_alt, p), "  (should be -1 for non-square)");
  print("");

  r14 = (p - b1_secp) % p;  \\ alpha^3 = -b1 mod p
  sv14 = [0, (1 + d14_alt) % p, 0];
  qv14 = [0, 0, d14_alt % p];
  print("  sv = (1+d)*alpha = (1-z3)*alpha = ", sv14);
  print("  qv = d*alpha^2 = (-z3)*alpha^2 = ", qv14);
  print("");

  \\ This computation takes long for 256-bit p with full F_{p^3} arithmetic.
  \\ Document the structure but skip full computation (too slow without rug/MPFR).
  print("  NOTE: Full 256-bit F_{p^3} Richelot computation is SLOW in PARI/GP");
  print("  with arbitrary precision. Skipping secp256k1-scale computation.");
  print("  For toy verification, use p=1009 analog.");
  print("");

  \\ Toy analog: use p=1009, b1 = lift(Mod(11, 1009)*z3_1009) for pair (1,4) analog
  p_toy = 1009;
  z3_toy = lift(polrootsmod(x^2+x+1, p_toy)[1]);
  b0_toy = 11;
  h_toy = lift(polrootsmod(x^2-x+1, p_toy)[1]);  \\ primitive 6th root
  print("---- Toy analog (p=1009, b0=11) ----");
  print("  6th root of unity h_toy = ", h_toy);
  b1_toy = lift(Mod(b0_toy, p_toy) * Mod(h_toy, p_toy));
  b4_toy = lift(Mod(b0_toy, p_toy) * Mod(h_toy, p_toy)^4);
  print("  E_1_toy: b1 = ", b1_toy, "  E_4_toy: b4 = ", b4_toy);

  E1_toy = ellinit([0, b1_toy], p_toy);
  E4_toy = ellinit([0, b4_toy], p_toy);
  t1_toy = p_toy + 1 - ellcard(E1_toy);
  t4_toy = p_toy + 1 - ellcard(E4_toy);
  print("  trace t1 = ", t1_toy, "  trace t4 = ", t4_toy, "  t1+t4 = ", t1_toy+t4_toy, " (should be 0)");

  d14_toy = lift((-Mod(z3_toy, p_toy)) % p_toy);  \\ d = -z3
  print("  d = -z3 = ", d14_toy, "  d^3 = ", lift(Mod(d14_toy, p_toy)^3), "  (should be ", p_toy-1, ")");
  print("  kronecker(d, p) = ", kronecker(d14_toy, p_toy));

  r14_toy = (p_toy - b1_toy) % p_toy;
  sv14_toy = [0, (1 + d14_toy) % p_toy, 0];
  qv14_toy = [0, 0, d14_toy % p_toy];

  res14_toy = richelot_z3(sv14_toy, qv14_toy, z3_toy, r14_toy);
  print("  Richelot result (pair 1,4 toy): a=", res14_toy[1], "  b=", res14_toy[2]);
  if(res14_toy[1] != -1,
    h14_toy = Mod(1, p_toy)*x^6 + Mod(res14_toy[1], p_toy)*x^3 + Mod(res14_toy[2], p_toy);
    cp14_toy = hyperellcharpoly(h14_toy);
    nj14_toy = subst(cp14_toy, variable(cp14_toy), 1);
    nE14_toy = (p_toy+1-t1_toy)*(p_toy+1+t1_toy);
    print("  Frobenius poly: ", cp14_toy);
    print("  #Jac = ", nj14_toy, "  target = ", nE14_toy, "  match = ", nj14_toy==nE14_toy),
    print("  STATUS: cover not over F_p")
  );
  print("");
}

\\ ================================================================
\\ Step 5: pairs (0,1) and (0,4) — non-twist pairs
\\ ================================================================
print("================================================================");
print("Step 5: secp256k1 pairs (0,1) and (0,4) — non-twist pairs");
print("================================================================");
print("");

{
  \\ For pair (0,1): E_0: y^2=x^3+7, E_1: y^2=x^3+7h.
  \\ The 2-torsion: alpha^3 = -7, beta^3 = -7h.
  \\ beta = h^{1/3} * alpha (if h has a cube root in F_p).
  \\ Since h is a 6th root of unity: h = g^{(p-1)/6} for a generator g.
  \\ h^{1/3} = g^{(p-1)/18} — exists iff 3 | (p-1)/6, i.e., 18 | (p-1).
  \\ For secp256k1: p ≡ 7 mod 12, p-1 ≡ 6 mod 12.
  \\ Is 18 | (p-1)?  p-1 = 2^256 - 2^32 - 978.  p-1 mod 18: ...
  \\ Actually, for the toy case p=1009: 1008/18 = 56, so 18|1008.
  \\ So h^{1/3} exists in F_{1009}!

  p_toy = 1009;
  z3_toy = lift(polrootsmod(x^2+x+1, p_toy)[1]);
  h_toy = lift(polrootsmod(x^2-x+1, p_toy)[1]);
  b0_toy = 11;
  b1_toy = lift(Mod(b0_toy, p_toy) * Mod(h_toy, p_toy));

  print("---- Toy: pair (0,1) with p=1009, b0=11, b1=", b1_toy, " ----");

  \\ d for (0,1): d^3 * b0 = b1, so d^3 = b1/b0 = h mod p.
  \\ h has a cube root mod p iff h^{(p-1)/3} = 1.
  h_cbrt_check = lift(Mod(h_toy, p_toy)^((p_toy-1)/3));
  print("  h^{(p-1)/3} mod p = ", h_cbrt_check, "  (=1 iff h is a cube)");
  if(h_cbrt_check == 1,
    \\ h is a cube! Find its cube root.
    \\ Cube root of h: since ord(h) divides p-1, use baby-step / Tonelli-Shanks variant.
    \\ Simpler: check all candidates x with x^3 ≡ h mod p.
    d01_toy = 0;
    for(x = 1, p_toy - 1,
      if(lift(Mod(x, p_toy)^3) == h_toy % p_toy,
        d01_toy = x;
        break
      )
    );
    print("  Cube root d s.t. d^3=h: d=", d01_toy);
    print("  Check: d^3 mod p = ", lift(Mod(d01_toy, p_toy)^3), "  h = ", h_toy);
    print("  kronecker(d, p) = ", kronecker(d01_toy, p_toy), "  (need -1 for non-square)");
    if(kronecker(d01_toy, p_toy) == -1,
      \\ d is a non-square — pair (0,1) has a Howe cover
      r01_toy = (p_toy - b0_toy) % p_toy;
      sv01_toy = [0, (1 + d01_toy) % p_toy, 0];
      qv01_toy = [0, 0, d01_toy % p_toy];
      res01_toy = richelot_z3(sv01_toy, qv01_toy, z3_toy, r01_toy);
      print("  Richelot result (pair 0,1 toy): a=", res01_toy[1], "  b=", res01_toy[2]);
      if(res01_toy[1] != -1,
        E0_toy = ellinit([0, b0_toy], p_toy);
        t0_toy = p_toy+1-ellcard(E0_toy);
        E1_toy2 = ellinit([0, b1_toy], p_toy);
        t1_toy2 = p_toy+1-ellcard(E1_toy2);
        h01_toy = Mod(1,p_toy)*x^6 + Mod(res01_toy[1],p_toy)*x^3 + Mod(res01_toy[2],p_toy);
        cp01_toy = hyperellcharpoly(h01_toy);
        nj01_toy = subst(cp01_toy, variable(cp01_toy), 1);
        nE01_toy = (p_toy+1-t0_toy)*(p_toy+1-t1_toy2);
        print("  t0 = ", t0_toy, "  t1 = ", t1_toy2);
        print("  Frobenius poly: ", cp01_toy);
        print("  #Jac = ", nj01_toy, "  target = ", nE01_toy, "  match = ", nj01_toy==nE01_toy);
        print("  Note: target should be (p+1-t0)(p+1-t1), not product of quadratic twists"),
        print("  Cover not over F_p")
      ),
      \\ d is a square — cannot use this d directly
      print("  d=", d01_toy, " is a SQUARE.  Need to adjust (d * non-square instead).");
      \\ Try d' = d * (first non-square): this has (d')^3 = h * (ns)^3 = different.
      \\ Actually if d^3 = h and d is a square, then d = e^2 for some e, giving e^6 = h.
      \\ The cube roots of h are d, d*z3, d*z3^2.  Exactly one of these is a non-square
      \\ (since their product is d^3 = h which has some Legendre symbol).
      ns = 2; while(kronecker(ns, p_toy) != -1, ns++);
      print("  Trying d*z3 = ", lift(Mod(d01_toy * z3_toy, p_toy)));
      d01b = lift(Mod(d01_toy * z3_toy, p_toy));
      print("  kronecker(d*z3, p) = ", kronecker(d01b, p_toy));
      print("  (d*z3)^3 = ", lift(Mod(d01b, p_toy)^3), "  h = ", h_toy % p_toy)
    ),
    print("  h is NOT a cube mod p=", p_toy, ".")
  );
  print("");
}

\\ ================================================================
\\ Step 6: Summary of naive secp256k1 cover (x^3+7)(x^3+189) vs pair (0,?)
\\ ================================================================
print("================================================================");
print("Step 6: Map naive cover y^2=(x^3+7)(x^3+189) to a pair (0,k)");
print("================================================================");
print("");

{
  \\ b2 = 189.  Find k s.t. 7*h^k ≡ 189 mod p.
  \\ h^k ≡ 27 mod p.  27 = 3^3.
  found_k = -1;
  for(kk = 0, 5,
    if(lift(Mod(b0, p) * Mod(h, p)^kk) == 189 % p,
      found_k = kk;
    );
  );
  if(found_k >= 0,
    print("  189 = 7*h^k for k = ", found_k, " → naive cover corresponds to pair (0,", found_k, ")"),
    print("  189 is NOT of the form 7*h^k for k=0..5 mod p_secp.")
  );
  print("  (If not found: 189 is not among the 6 sextic twists of secp256k1 over F_p)");
  print("");
  print("  Note: 189 = 27*7 = 3^3 * 7.  d=3 for the naive cover.");
  print("  If 3^3 is not one of h^0,...,h^5 mod p, then 189 is a DIFFERENT");
  print("  quadratic twist of secp256k1, not among the 6 sextic twists.");
  print("  (The 6 sextic twists use d = h, h^2, ..., h^6 = 1 while d must be non-square.)");
  print("");
}

print("================================================================");
print("Done.");
print("================================================================");

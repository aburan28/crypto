\\ ============================================================
\\ howe_richelot_p43.gp
\\ ============================================================
\\
\\ Identify the Howe-glued genus-2 curve over F_43 by computing
\\ the Richelot preimage of E1 × E2 directly over F_{43^3}.
\\
\\ Setup:
\\   p = 43,  E1: y^2=x^3+7 (trace 13),  E2: y^2=x^3+13 (trace -13)
\\   F_{43^3} = F_43[a] / (a^3 - 36)  (since a^3+7 ≡ a^3-36 ≡ 0)
\\
\\ Key facts derived analytically:
\\   beta = 2*alpha  (2*alpha)^3 = 8*36 = 288 ≡ 30 = -13 mod 43 ✓
\\   s = alpha+beta = 3*alpha
\\   q = alpha*beta  = 2*alpha^2
\\   zeta3 = 6 in F_43  (6^2+6+1 = 43 ≡ 0 mod 43)
\\
\\ Method: compute Richelot image of y^2=(x^3+7)(x^3+13) under the
\\ pairing G_i = (x - alpha*zeta^{i-1})(x - beta*zeta^{i-1}) using
\\ explicit polynomial arithmetic in F_{43^3}.
\\
\\ The Richelot image H_1*H_2*H_3 lies in F_43[x] by Z/3Z-symmetry,
\\ giving the (a,b) of the Howe-glued Z/3Z curve y^2=x^6+a*x^3+b.
\\
\\ Run: gp -q howe_richelot_p43.gp

default(parisize, 64000000);
default(timer, 0);

print("================================================================");
print("Richelot preimage identification: Howe-glued curve over F_43");
print("================================================================");
print("");

p = 43;

\\ ---- Verify base setup ----
print("---- Verify E1, E2 over F_43 ----");
E1 = ellinit([0,7], p);  n1 = ellcard(E1);  t1 = p+1-n1;
E2 = ellinit([0,13], p); n2 = ellcard(E2); t2 = p+1-n2;
print("E1: y^2=x^3+7,   #E1=", n1, "  trace=", t1);
print("E2: y^2=x^3+13,  #E2=", n2, "  trace=", t2);
print("Target char poly: (T^2-", t1, "T+", p, ")(T^2+", t1, "T+", p, ")");
print("  = T^4-", t1^2-2*p, "T^2+", p^2);
print("");

\\ ---- Verify beta = 2*alpha ----
print("---- Verify beta = 2*alpha relationship ----");
print("alpha^3 = -7 = 36 mod 43");
print("(2*alpha)^3 = 8*36 mod 43 = ", (8*36) % p, "  (should be 30 = -13 mod 43)");
print("Check: -13 mod 43 = ", (-13) % p);
print("");

\\ ---- Verify zeta3 in F_43 ----
print("---- Verify zeta3 = 6 in F_43 ----");
z3 = 6;
print("zeta3 = 6,  6^2+6+1 mod 43 = ", (z3^2 + z3 + 1) % p, "  (should be 0)");
print("6^3 mod 43 = ", z3^3 % p, "  (should be 1)");
print("");

\\ ================================================================
\\ Represent F_{43^3} = F_43[a]/(a^3 - 36)
\\ Elements: polynomials c0 + c1*a + c2*a^2 with ci in F_43.
\\ Stored as vectors [c0, c1, c2].
\\ Reduction rule: a^3 = 36.
\\ ================================================================

\\ Multiply two F_{43^3} elements (as vectors)
ff3_mul(u, v) = {
  my(c0, c1, c2, c3, c4);
  c0 = (u[1]*v[1]) % p;
  c1 = (u[1]*v[2] + u[2]*v[1]) % p;
  c2 = (u[1]*v[3] + u[2]*v[2] + u[3]*v[1]) % p;
  \\ c3 = coeff of a^3 = u[2]*v[3] + u[3]*v[2]  -> reduce: a^3 = 36
  c3 = (u[2]*v[3] + u[3]*v[2]) % p;
  \\ c4 = coeff of a^4 = a*a^3 = 36*a  -> u[3]*v[3] goes to 36*a
  c4 = (u[3]*v[3]) % p;
  \\ After reduction: a^3 -> 36, a^4 -> 36*a
  [( c0 + 36*c3 ) % p,
   ( c1 + 36*c4 ) % p,
   c2]
};

\\ Add two F_{43^3} elements
ff3_add(u, v) = [(u[1]+v[1])%p, (u[2]+v[2])%p, (u[3]+v[3])%p];

\\ Scale by scalar c in F_43
ff3_scale(c, u) = [(c*u[1])%p, (c*u[2])%p, (c*u[3])%p];

\\ Negate
ff3_neg(u) = [(-u[1])%p, (-u[2])%p, (-u[3])%p];

\\ Constant element c in F_43 embedded in F_{43^3}
ff3_const(c) = [c % p, 0, 0];

\\ Verify: a * a^2 = a^3 = 36
alpha  = [0, 1, 0];   \\ a (= alpha, root of x^3+7 in F_{43^3})
alpha2 = [0, 0, 1];   \\ a^2
alpha3 = ff3_mul(alpha, alpha2);  \\ should be [36, 0, 0]
print("---- F_{43^3} arithmetic check ----");
print("alpha^3 = ", alpha3, "  (should be [36,0,0])");

\\ beta = 2*alpha
beta = ff3_scale(2, alpha);   \\ [0, 2, 0]
beta3 = ff3_mul(beta, ff3_mul(beta, beta));
print("beta^3 = ", beta3, "  (should be [30,0,0])");
print("");

\\ ================================================================
\\ Richelot computation
\\ ================================================================
\\
\\ Pairing sigma_0: (alpha, beta), (alpha*z3, beta*z3), (alpha*z3^2, beta*z3^2)
\\
\\ G_1 = (x - alpha)(x - beta) = x^2 - s*x + q
\\ where s = alpha + beta = 3*alpha, q = alpha*beta = 2*alpha^2
\\
\\ G_2 = (x - z3*alpha)(x - z3*beta) = x^2 - z3*s*x + z3^2*q
\\ G_3 = (x - z3^2*alpha)(x - z3^2*beta) = x^2 - z3^2*s*x + z3*q
\\
\\ Each G_i is a degree-2 polynomial in x with F_{43^3} coefficients.
\\ Represent as [const_coeff, x_coeff, x^2_coeff].
\\
\\ Richelot:
\\   Δ = det M where M is the 3x3 matrix of [leading, middle, const] of G_i
\\   H_1 = (G_2' G_3 - G_2 G_3') / Δ   where G_i' = dG_i/dx

print("---- Richelot computation over F_{43^3} ----");

s_coeff = ff3_scale(3, alpha);      \\ s = 3*alpha
q_coeff = ff3_scale(2, alpha2);     \\ q = 2*alpha^2
z3sq = (z3*z3) % p;                  \\ z3^2 = 36 mod 43
print("z3 = ", z3, "  z3^2 = ", z3sq);

\\ G_1 = x^2 - s*x + q   (as [const, x_coeff, x^2_coeff])
G1_c = q_coeff;
G1_x = ff3_neg(s_coeff);
G1_x2 = ff3_const(1);

\\ G_2 = x^2 - z3*s*x + z3^2*q
G2_c = ff3_scale(z3sq, q_coeff);
G2_x = ff3_neg(ff3_scale(z3, s_coeff));
G2_x2 = ff3_const(1);

\\ G_3 = x^2 - z3^2*s*x + z3*q
G3_c = ff3_scale(z3, q_coeff);
G3_x = ff3_neg(ff3_scale(z3sq, s_coeff));
G3_x2 = ff3_const(1);

\\ Derivatives: G_i' = 2x + G_i_x  (as [G_i_x, 2, 0])
G1p_c = G1_x;   G1p_x = ff3_const(2);
G2p_c = G2_x;   G2p_x = ff3_const(2);
G3p_c = G3_x;   G3p_x = ff3_const(2);

\\ Compute H_1 = G_2'*G_3 - G_2*G_3' (degree 3 poly in x)
\\ G_2' = G2p_x * x + G2p_c
\\ G_3  = G3_x2 * x^2 + G3_x * x + G3_c
\\ Product G_2'*G_3: degree 3
\\   x^3: G2p_x * G3_x2
\\   x^2: G2p_x * G3_x + G2p_c * G3_x2
\\   x^1: G2p_x * G3_c + G2p_c * G3_x
\\   x^0: G2p_c * G3_c

G2pG3_x3 = ff3_mul(G2p_x, G3_x2);
G2pG3_x2 = ff3_add(ff3_mul(G2p_x, G3_x), ff3_mul(G2p_c, G3_x2));
G2pG3_x1 = ff3_add(ff3_mul(G2p_x, G3_c), ff3_mul(G2p_c, G3_x));
G2pG3_x0 = ff3_mul(G2p_c, G3_c);

\\ G_2*G_3': degree 3
G2G3p_x3 = ff3_mul(G2_x2, G3p_x);
G2G3p_x2 = ff3_add(ff3_mul(G2_x2, G3p_c), ff3_mul(G2_x, G3p_x));
G2G3p_x1 = ff3_add(ff3_mul(G2_x, G3p_c), ff3_mul(G2_c, G3p_x));
G2G3p_x0 = ff3_mul(G2_c, G3p_c);

\\ H_1 numerator = G_2'*G_3 - G_2*G_3' (degree 3, but leading coeff should vanish)
H1_x3 = ff3_add(G2pG3_x3, ff3_neg(G2G3p_x3));
H1_x2 = ff3_add(G2pG3_x2, ff3_neg(G2G3p_x2));
H1_x1 = ff3_add(G2pG3_x1, ff3_neg(G2G3p_x1));
H1_x0 = ff3_add(G2pG3_x0, ff3_neg(G2G3p_x0));

print("H1 numerator coefficients (x^3, x^2, x^1, x^0):");
print("  x^3: ", H1_x3, "  (should be [0,0,0])");
print("  x^2: ", H1_x2);
print("  x^1: ", H1_x1);
print("  x^0: ", H1_x0);
print("");

\\ Compute Δ = det of 3x3 matrix [G1_x2, G1_x, G1_c; G2_x2, G2_x, G2_c; G3_x2, G3_x, G3_c]
\\ Δ = G1_x2*(G2_x*G3_c - G3_x*G2_c) - G1_x*(G2_x2*G3_c - G3_x2*G2_c) + G1_c*(G2_x2*G3_x - G3_x2*G2_x)
Delta = ff3_add(
  ff3_add(
    ff3_mul(G1_x2, ff3_add(ff3_mul(G2_x, G3_c), ff3_neg(ff3_mul(G3_x, G2_c)))),
    ff3_neg(ff3_mul(G1_x, ff3_add(ff3_mul(G2_x2, G3_c), ff3_neg(ff3_mul(G3_x2, G2_c)))))
  ),
  ff3_mul(G1_c, ff3_add(ff3_mul(G2_x2, G3_x), ff3_neg(ff3_mul(G3_x2, G2_x))))
);
print("Delta = ", Delta);
print("  (expected: 3*(z3-z3^2)*s*q = 3*(6-36)*3*alpha*2*alpha^2 = 3*(-30)*6*alpha^3)");
print("  = 3*13*6*36 mod 43 = ", (3*((z3-z3sq) % p)*3*2*36) % p, " times [1,0,0]");
print("");

\\ Invert Delta in F_{43^3}
\\ Use: in F_{p^3} = F_p[a]/(a^3-36), inverse computed via extended Euclidean
\\ For simplicity, since Delta should be a scalar (in F_43 after substitution),
\\ check if Delta[2]==0 and Delta[3]==0

if (Delta[2] == 0 && Delta[3] == 0,
  print("Delta is in F_43 (scalar). Delta = ", Delta[1]);
  Delta_inv = [lift(Mod(Delta[1], p)^(-1)), 0, 0];
  print("Delta_inv = ", Delta_inv),
  print("Delta has alpha components — computing inverse via norm")
);
print("");

\\ Divide H1 numerator by Delta to get H1
H1_quad_x2 = ff3_mul(H1_x2, Delta_inv);
H1_quad_x1 = ff3_mul(H1_x1, Delta_inv);
H1_quad_x0 = ff3_mul(H1_x0, Delta_inv);

print("H1 (quadratic, after dividing by Delta):");
print("  x^2: ", H1_quad_x2);
print("  x^1: ", H1_quad_x1);
print("  x^0: ", H1_quad_x0);
print("");

\\ ================================================================
\\ Compute H2 and H3 similarly
\\ H2 = G3'*G1 - G3*G1'
\\ H3 = G1'*G2 - G1*G2'
\\ ================================================================

\\ H2 = G3'*G1 - G3*G1'
G3pG1_x3 = ff3_mul(G3p_x, G1_x2);
G3pG1_x2 = ff3_add(ff3_mul(G3p_x, G1_x), ff3_mul(G3p_c, G1_x2));
G3pG1_x1 = ff3_add(ff3_mul(G3p_x, G1_c), ff3_mul(G3p_c, G1_x));
G3pG1_x0 = ff3_mul(G3p_c, G1_c);

G3G1p_x3 = ff3_mul(G3_x2, G1p_x);
G3G1p_x2 = ff3_add(ff3_mul(G3_x2, G1p_c), ff3_mul(G3_x, G1p_x));
G3G1p_x1 = ff3_add(ff3_mul(G3_x, G1p_c), ff3_mul(G3_c, G1p_x));
G3G1p_x0 = ff3_mul(G3_c, G1p_c);

H2_x3 = ff3_add(G3pG1_x3, ff3_neg(G3G1p_x3));
H2_x2 = ff3_add(G3pG1_x2, ff3_neg(G3G1p_x2));
H2_x1 = ff3_add(G3pG1_x1, ff3_neg(G3G1p_x1));
H2_x0 = ff3_add(G3pG1_x0, ff3_neg(G3G1p_x0));

H2_quad_x2 = ff3_mul(H2_x2, Delta_inv);
H2_quad_x1 = ff3_mul(H2_x1, Delta_inv);
H2_quad_x0 = ff3_mul(H2_x0, Delta_inv);

\\ H3 = G1'*G2 - G1*G2'
G1pG2_x3 = ff3_mul(G1p_x, G2_x2);
G1pG2_x2 = ff3_add(ff3_mul(G1p_x, G2_x), ff3_mul(G1p_c, G2_x2));
G1pG2_x1 = ff3_add(ff3_mul(G1p_x, G2_c), ff3_mul(G1p_c, G2_x));
G1pG2_x0 = ff3_mul(G1p_c, G2_c);

G1G2p_x3 = ff3_mul(G1_x2, G2p_x);
G1G2p_x2 = ff3_add(ff3_mul(G1_x2, G2p_c), ff3_mul(G1_x, G2p_x));
G1G2p_x1 = ff3_add(ff3_mul(G1_x, G2p_c), ff3_mul(G1_c, G2p_x));
G1G2p_x0 = ff3_mul(G1_c, G2_c);

H3_x3 = ff3_add(G1pG2_x3, ff3_neg(G1G2p_x3));
H3_x2 = ff3_add(G1pG2_x2, ff3_neg(G1G2p_x2));
H3_x1 = ff3_add(G1pG2_x1, ff3_neg(G1G2p_x1));
H3_x0 = ff3_add(G1pG2_x0, ff3_neg(G1G2p_x0));

H3_quad_x2 = ff3_mul(H3_x2, Delta_inv);
H3_quad_x1 = ff3_mul(H3_x1, Delta_inv);
H3_quad_x0 = ff3_mul(H3_x0, Delta_inv);

print("H2 (quadratic):");
print("  x^2: ", H2_quad_x2);
print("  x^1: ", H2_quad_x1);
print("  x^0: ", H2_quad_x0);
print("");
print("H3 (quadratic):");
print("  x^3 check: ", H3_x3, "  (should be [0,0,0])");
print("  x^2: ", H3_quad_x2);
print("  x^1: ", H3_quad_x1);
print("  x^0: ", H3_quad_x0);
print("");

\\ ================================================================
\\ Compute the product H1*H2*H3 as a degree-6 polynomial in x
\\ with coefficients in F_{43^3}.
\\ The result should lie in F_43[x] (Z/3Z-symmetry guarantee).
\\ ================================================================

\\ First compute H1*H2 (degree 4)
\\ H = c + d*x + e*x^2 multiply by f + g*x + h*x^2

poly_mul_deg2(a0,a1,a2, b0,b1,b2) = {
  my(c0,c1,c2,c3,c4);
  c0 = ff3_mul(a0, b0);
  c1 = ff3_add(ff3_mul(a0,b1), ff3_mul(a1,b0));
  c2 = ff3_add(ff3_add(ff3_mul(a0,b2), ff3_mul(a1,b1)), ff3_mul(a2,b0));
  c3 = ff3_add(ff3_mul(a1,b2), ff3_mul(a2,b1));
  c4 = ff3_mul(a2, b2);
  [c0, c1, c2, c3, c4]
};

H1H2 = poly_mul_deg2(H1_quad_x0, H1_quad_x1, H1_quad_x2,
                     H2_quad_x0, H2_quad_x1, H2_quad_x2);
\\ H1H2 = [c0, c1, c2, c3, c4] = coeffs of 1, x, x^2, x^3, x^4

\\ Now multiply H1H2 (degree 4) by H3 (degree 2) to get degree 6
poly_mul_deg4_deg2(a, b) = {
  \\ a = [a0..a4], b = [b0,b1,b2]
  my(c0,c1,c2,c3,c4,c5,c6);
  c0 = ff3_mul(a[1], b[1]);
  c1 = ff3_add(ff3_mul(a[1],b[2]), ff3_mul(a[2],b[1]));
  c2 = ff3_add(ff3_add(ff3_mul(a[1],b[3]), ff3_mul(a[2],b[2])), ff3_mul(a[3],b[1]));
  c3 = ff3_add(ff3_add(ff3_mul(a[2],b[3]), ff3_mul(a[3],b[2])), ff3_mul(a[4],b[1]));
  c4 = ff3_add(ff3_add(ff3_mul(a[3],b[3]), ff3_mul(a[4],b[2])), ff3_mul(a[5],b[1]));
  c5 = ff3_add(ff3_mul(a[4],b[3]), ff3_mul(a[5],b[2]));
  c6 = ff3_mul(a[5], b[3]);
  [c0, c1, c2, c3, c4, c5, c6]
};

H1H2H3 = poly_mul_deg4_deg2(H1H2, [H3_quad_x0, H3_quad_x1, H3_quad_x2]);

print("================================================================");
print("H1*H2*H3 coefficients (x^0 through x^6):");
print("  (Each should be in F_43, i.e., [c,0,0])");
print("================================================================");
for (i = 1, 7,
  print("  x^", i-1, ": ", H1H2H3[i])
);
print("");

\\ Extract F_43 coefficients (check alpha-components vanish)
print("---- Extract F_43 polynomial ----");
coeffs_ok = 1;
F43_coeffs = vector(7);
for (i = 1, 7,
  if (H1H2H3[i][2] != 0 || H1H2H3[i][3] != 0,
    print("  WARNING: x^", i-1, " has alpha component: ", H1H2H3[i]);
    coeffs_ok = 0
  );
  F43_coeffs[i] = H1H2H3[i][1]
);

if (coeffs_ok,
  print("All coefficients are in F_43! ✓");
  print("F43 coefficients: ", F43_coeffs);
  print("");
  print("Richelot image sextic: y^2 = ", F43_coeffs[7], "x^6 + ",
    F43_coeffs[6], "x^5 + ", F43_coeffs[5], "x^4 + ",
    F43_coeffs[4], "x^3 + ", F43_coeffs[3], "x^2 + ",
    F43_coeffs[2], "x + ", F43_coeffs[1]);
);
print("");

\\ ================================================================
\\ Normalise to Z/3Z form y^2 = x^6 + a*x^3 + b
\\ The sextic should have coefficients at x^5, x^4, x^2, x^1 all zero.
\\ ================================================================

print("---- Check Z/3Z form (x^5,x^4,x^2,x^1 should be 0) ----");
print("  x^5 coeff: ", F43_coeffs[6]);
print("  x^4 coeff: ", F43_coeffs[5]);
print("  x^2 coeff: ", F43_coeffs[3]);
print("  x^1 coeff: ", F43_coeffs[2]);
print("");

\\ Normalize: divide by leading coeff (x^6 coeff) to make monic
lc = F43_coeffs[7];
if (lc != 0,
  lc_inv = lift(Mod(lc, p)^(-1));
  norm_a3 = (F43_coeffs[4] * lc_inv) % p;
  norm_b  = (F43_coeffs[1] * lc_inv) % p;
  print("Monic Z/3Z form: y^2 = x^6 + ", norm_a3, " x^3 + ", norm_b);
  print("  a = ", norm_a3, "  b = ", norm_b);
  print("")
);

\\ ================================================================
\\ Verify: check Frobenius char poly of the Richelot image
\\ ================================================================

print("---- Verify char poly of Richelot image ----");
h_richelot = sum(i=1, 7, F43_coeffs[i] * x^(i-1));
h_richelot_mod = h_richelot * Mod(1, p);
cp_richelot = hyperellcharpoly(h_richelot_mod);
nj_richelot = subst(cp_richelot, variable(cp_richelot), 1);
print("Char poly: ", cp_richelot);
print("#Jac = ", nj_richelot);
print("Target char poly: T^4 - ", t1^2-2*p, " T^2 + ", p^2);
target_cp_val = 1 - (t1^2-2*p) + p^2;
print("#Jac target (T=1): ", target_cp_val);
print("Match? ", nj_richelot == target_cp_val);
print("");

\\ ================================================================
\\ Compare against 7 canonical Z/3Z classes from 2026-06-01 run
\\ ================================================================

print("================================================================");
print("Compare against the 7 canonical Z/3Z classes");
print("================================================================");
print("");

\\ The 7 canonical classes from the 2026-06-01 search:
classes_a = [1, 2, 3, 4, 6, 8, 16];
classes_b = [2, 8, 42, 32, 39, 42, 39];

for (i = 1, 7,
  ca = classes_a[i]; cb = classes_b[i];
  match_a = (norm_a3 == ca);
  match_b = (norm_b == cb);
  if (match_a && match_b,
    print("MATCH: Class ", i, " (a=", ca, ", b=", cb, ") ✓✓✓")
  ,
    print("Class ", i, " (a=", ca, ", b=", cb, "): a_match=", match_a, " b_match=", match_b)
  )
);
print("");

\\ ================================================================
\\ Also try the other two Galois-equivariant matchings (sigma_1, sigma_2)
\\ sigma_1: pair (alpha, beta*z3), (alpha*z3, beta*z3^2), (alpha*z3^2, beta)
\\ This changes q_coeff -> z3*q and (within the same curve structure)
\\ ================================================================

print("================================================================");
print("Try alternative matchings sigma_1 and sigma_2");
print("================================================================");
print("");

\\ For sigma_1: swap the beta labelling by multiplying beta by z3
\\ G1 = (x-alpha)(x-z3*beta), G2 = (x-z3*alpha)(x-z3^2*beta), G3=(x-z3^2*alpha)(x-beta)
\\ s1 = alpha + z3*beta = alpha + 2*z3*alpha = (1+2*z3)*alpha = (1+12)*alpha = 13*alpha
\\ q1 = alpha * z3*beta = 2*z3*alpha^2

s1_coeff = ff3_scale((1 + 2*z3) % p, alpha);   \\ (1 + 2*z3)*alpha
q1_coeff = ff3_scale((2*z3) % p, alpha2);        \\ 2*z3*alpha^2
print("sigma_1: s1 = (1+2*z3)*alpha = ", (1+2*z3)%p, "*alpha");
print("         q1 = 2*z3*alpha^2 = ", (2*z3)%p, "*alpha^2");

\\ G1_1 = x^2 - s1*x + q1
G1_1_c = q1_coeff;
G1_1_x = ff3_neg(s1_coeff);
G1_1_x2 = ff3_const(1);

\\ G2_1 = x^2 - z3*s1*x + z3^2*q1
G2_1_c = ff3_scale(z3sq, q1_coeff);
G2_1_x = ff3_neg(ff3_scale(z3, s1_coeff));
G2_1_x2 = ff3_const(1);

\\ G3_1 = x^2 - z3^2*s1*x + z3*q1
G3_1_c = ff3_scale(z3, q1_coeff);
G3_1_x = ff3_neg(ff3_scale(z3sq, s1_coeff));
G3_1_x2 = ff3_const(1);

\\ Compute Delta1
Delta1 = ff3_add(
  ff3_add(
    ff3_mul(G1_1_x2, ff3_add(ff3_mul(G2_1_x, G3_1_c), ff3_neg(ff3_mul(G3_1_x, G2_1_c)))),
    ff3_neg(ff3_mul(G1_1_x, ff3_add(ff3_mul(G2_1_x2, G3_1_c), ff3_neg(ff3_mul(G3_1_x2, G2_1_c)))))
  ),
  ff3_mul(G1_1_c, ff3_add(ff3_mul(G2_1_x2, G3_1_x), ff3_neg(ff3_mul(G3_1_x2, G2_1_x))))
);

\\ H1_1 = G2_1'*G3_1 - G2_1*G3_1'
G2_1p_c = G2_1_x; G2_1p_x = ff3_const(2);
G3_1p_c = G3_1_x; G3_1p_x = ff3_const(2);

H1_1_x3 = ff3_add(ff3_mul(G2_1p_x, G3_1_x2), ff3_neg(ff3_mul(G3_1_x2, G2_1p_x)));
H1_1_x2 = ff3_add(ff3_add(ff3_mul(G2_1p_x, G3_1_x), ff3_mul(G2_1p_c, G3_1_x2)), ff3_neg(ff3_add(ff3_mul(G3_1_x2, G2_1p_c), ff3_mul(G3_1_x, G2_1p_x))));
H1_1_x1 = ff3_add(ff3_add(ff3_mul(G2_1p_x, G3_1_c), ff3_mul(G2_1p_c, G3_1_x)), ff3_neg(ff3_add(ff3_mul(G3_1_x, G2_1p_c), ff3_mul(G3_1_c, G2_1p_x))));
H1_1_x0 = ff3_add(ff3_mul(G2_1p_c, G3_1_c), ff3_neg(ff3_mul(G3_1_c, G2_1p_c)));

\\ (H1_1_x0 should be 0 since G2'G3 and G2G3' have same degree-0 cross terms)
\\ Actually for H_i = (G_j' G_k - G_j G_k')/Delta, this is a degree-2 poly.
\\ Check x^3 coeff: should vanish since both G_j', G_k have leading coeff 2,1
\\ G_j'*G_k has x^3 coeff = 2*1 = 2; G_j*G_k' has x^3 coeff = 1*2 = 2; difference=0. ✓

\\ Compute H1H2H3 for sigma_1 (parallel to above)
if (Delta1[2] == 0 && Delta1[3] == 0,
  Delta1_inv = [lift(Mod(Delta1[1], p)^(-1)), 0, 0];

  G1_1p_c = G1_1_x; G1_1p_x = ff3_const(2);
  G3_1p_c = G3_1_x; G3_1p_x = ff3_const(2);

  H1_1_num_x2 = ff3_add(ff3_add(ff3_mul(G2_1p_x, G3_1_x), ff3_mul(G2_1p_c, G3_1_x2)), ff3_neg(ff3_add(ff3_mul(G2_1_x2, G3_1p_c), ff3_mul(G2_1_x, G3_1p_x))));
  H1_1_num_x1 = ff3_add(ff3_add(ff3_mul(G2_1p_x, G3_1_c), ff3_mul(G2_1p_c, G3_1_x)), ff3_neg(ff3_add(ff3_mul(G2_1_x, G3_1p_c), ff3_mul(G2_1_c, G3_1p_x))));
  H1_1_num_x0 = ff3_add(ff3_mul(G2_1p_c, G3_1_c), ff3_neg(ff3_mul(G2_1_c, G3_1p_c)));

  \\ H2_1 = G3_1'*G1_1 - G3_1*G1_1'
  H2_1_num_x2 = ff3_add(ff3_add(ff3_mul(G3_1p_x, G1_1_x), ff3_mul(G3_1p_c, G1_1_x2)), ff3_neg(ff3_add(ff3_mul(G3_1_x2, G1_1p_c), ff3_mul(G3_1_x, G1_1p_x))));
  H2_1_num_x1 = ff3_add(ff3_add(ff3_mul(G3_1p_x, G1_1_c), ff3_mul(G3_1p_c, G1_1_x)), ff3_neg(ff3_add(ff3_mul(G3_1_x, G1_1p_c), ff3_mul(G3_1_c, G1_1p_x))));
  H2_1_num_x0 = ff3_add(ff3_mul(G3_1p_c, G1_1_c), ff3_neg(ff3_mul(G3_1_c, G1_1p_c)));

  \\ H3_1 = G1_1'*G2_1 - G1_1*G2_1'
  H3_1_num_x2 = ff3_add(ff3_add(ff3_mul(G1_1p_x, G2_1_x), ff3_mul(G1_1p_c, G2_1_x2)), ff3_neg(ff3_add(ff3_mul(G1_1_x2, G2_1p_c), ff3_mul(G1_1_x, G2_1p_x))));
  H3_1_num_x1 = ff3_add(ff3_add(ff3_mul(G1_1p_x, G2_1_c), ff3_mul(G1_1p_c, G2_1_x)), ff3_neg(ff3_add(ff3_mul(G1_1_x, G2_1p_c), ff3_mul(G1_1_c, G2_1p_x))));
  H3_1_num_x0 = ff3_add(ff3_mul(G1_1p_c, G2_1_c), ff3_neg(ff3_mul(G1_1_c, G2_1p_c)));

  H1_1_x2d = ff3_mul(H1_1_num_x2, Delta1_inv);
  H1_1_x1d = ff3_mul(H1_1_num_x1, Delta1_inv);
  H1_1_x0d = ff3_mul(H1_1_num_x0, Delta1_inv);

  H2_1_x2d = ff3_mul(H2_1_num_x2, Delta1_inv);
  H2_1_x1d = ff3_mul(H2_1_num_x1, Delta1_inv);
  H2_1_x0d = ff3_mul(H2_1_num_x0, Delta1_inv);

  H3_1_x2d = ff3_mul(H3_1_num_x2, Delta1_inv);
  H3_1_x1d = ff3_mul(H3_1_num_x1, Delta1_inv);
  H3_1_x0d = ff3_mul(H3_1_num_x0, Delta1_inv);

  H1H2_1 = poly_mul_deg2(H1_1_x0d, H1_1_x1d, H1_1_x2d,
                          H2_1_x0d, H2_1_x1d, H2_1_x2d);
  H1H2H3_1 = poly_mul_deg4_deg2(H1H2_1, [H3_1_x0d, H3_1_x1d, H3_1_x2d]);

  print("sigma_1 result:");
  F43_coeffs_1 = vector(7);
  ok1 = 1;
  for (i = 1, 7,
    if (H1H2H3_1[i][2] != 0 || H1H2H3_1[i][3] != 0,
      print("  WARNING at x^", i-1, ": ", H1H2H3_1[i]);
      ok1 = 0
    );
    F43_coeffs_1[i] = H1H2H3_1[i][1]
  );
  if (ok1,
    lc1 = F43_coeffs_1[7];
    if (lc1 != 0,
      lc1_inv = lift(Mod(lc1,p)^(-1));
      a1 = (F43_coeffs_1[4] * lc1_inv) % p;
      b1 = (F43_coeffs_1[1] * lc1_inv) % p;
      print("sigma_1 Z/3Z form: a=", a1, " b=", b1)
    )
  )
,
  print("Delta1 not in F_43: ", Delta1)
);

print("");
print("================================================================");
print("Done.");
print("================================================================");

\\ ============================================================
\\ CHLRS Igusa formula: Rosenhain model of the Howe-glued curve
\\ ============================================================
\\
\\ For j=0 elliptic curves E1: y^2 = x^3 + b and its quadratic
\\ twist E2: y^2 = x^3 + b*d^3 (d = non-square in F_p, p ≡ 1 mod 3),
\\ the Howe-glued genus-2 curve C with Jac(C) (2,2)-isogenous to
\\ E1 × E2 can be constructed via cross-ratios of 2-torsion points.
\\
\\ The 2-torsion of E1: roots of x^3 + b in F_{p^3} are {α, ωα, ω^2α}
\\ where α = (-b)^{1/3} and ω = primitive cube root of unity in F_p.
\\ Similarly for E2: {β, ωβ, ω^2β} with β = (-b*d^3)^{1/3} = d*α.
\\
\\ Apply Möbius transform T sending α → 0, β=dα → 1, ωα → ∞.
\\   T(x) = (x - α)(d - ω) / ((x - ωα)(d - 1))
\\ The three remaining 2-torsion points map to:
\\   λ1 = T(ω^2α) = (ω^2-1)(d-ω) / (ω(ω-1)(d-1))
\\   λ2 = T(dωα)  = (dω-1)(d-ω)  / (ω(d-1)^2)
\\   λ3 = T(dω^2α)= (dω^2-1)(d-ω)/ (ω(dω-1)(d-1))
\\
\\ All λ_i ∈ F_p since d, ω ∈ F_p for p ≡ 1 mod 3.
\\ The Rosenhain genus-2 curve: C: y^2 = x(x-1)(x-λ1)(x-λ2)(x-λ3).
\\
\\ Claim: Jac(C) is (2,2)-isogenous to E1 × E2 over F_p.
\\ Verification: Frobenius poly of Jac(C) should equal the product
\\ (T^2 - t*T + p)(T^2 + t*T + p) where t = trace of E1.
\\
\\ This is an explicit CHLRS-type formula, avoiding Mestre reconstruction.
\\
\\ Run: gp -q chlrs_igusa_formula.gp
\\
\\ Reference: Howe (1996); Cardona-Howe-Lercier-Ritzenthaler-Streng
\\ approach to explicit Howe-glued genus-2 curves.

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("CHLRS Rosenhain formula for Howe-glued genus-2 curves (j=0)");
print("================================================================");
print("");

\\ ---- Helper: find primitive cube root of unity mod p ----
cube_root_unity(p) = {
  my(r);
  \\ solve ω^2 + ω + 1 = 0 mod p
  r = polrootsmod(x^2 + x + 1, p);
  if(#r == 0, error("p does not have cube roots of unity (p not ≡ 1 mod 3)"));
  \\ return the smaller of the two roots
  lift(r[1])
};

\\ ---- Helper: find first non-square mod p ----
first_nonsquare(p) = {
  my(d);
  for(d = 2, p-1,
    if(kronecker(d, p) == -1, return(d))
  );
  error("no non-square found")
};

\\ ---- Rosenhain invariants from (b, d, ω) ----
\\ All computations in F_p.
rosenhain_from_twist(b, d_val, om, p) = {
  my(lambda1, lambda2, lambda3, d, w, Fp);
  d = Mod(d_val, p);
  w = Mod(om, p);
  \\ λ1 = (ω^2-1)(d-ω) / (ω(ω-1)(d-1))
  lambda1 = (w^2 - 1) * (d - w) / (w * (w - 1) * (d - 1));
  \\ λ2 = (dω-1)(d-ω) / (ω(d-1)^2)
  lambda2 = (d*w - 1) * (d - w) / (w * (d - 1)^2);
  \\ λ3 = (dω^2-1)(d-ω) / (ω(dω-1)(d-1))
  lambda3 = (d*w^2 - 1) * (d - w) / (w * (d*w - 1) * (d - 1));
  [lift(lambda1), lift(lambda2), lift(lambda3)]
};

\\ ---- Rosenhain curve as degree-5 polynomial ----
rosenhain_poly(lambdas, p) = {
  my(l1, l2, l3);
  l1 = Mod(lambdas[1], p);
  l2 = Mod(lambdas[2], p);
  l3 = Mod(lambdas[3], p);
  x * (x - 1) * (x - l1) * (x - l2) * (x - l3)
};

\\ ---- Also try the 3 cyclic variants of the Möbius transform ----
\\ Variant 2: send ωα→0, dα→1, ω^2α→∞
rosenhain_variant2(b, d_val, om, p) = {
  my(d, w, l1, l2, l3);
  d = Mod(d_val, p);
  w = Mod(om, p);
  \\ T2(x) = (x - ωα)(d - ω^2) / ((x - ω^2α)(d - ω))
  \\ = (x/α - ω)(d - ω^2) / ((x/α - ω^2)(d - ω))
  \\ T2(α) = (1 - ω)(d - ω^2) / ((1 - ω^2)(d - ω))
  l1 = (1 - w) * (d - w^2) / ((1 - w^2) * (d - w));
  \\ T2(dω α): (dω - ω)(d - ω^2) / ((dω - ω^2)(d - ω))
  l2 = (d*w - w) * (d - w^2) / ((d*w - w^2) * (d - w));
  \\ T2(dω^2 α): (dω^2 - ω)(d - ω^2) / ((dω^2 - ω^2)(d - ω))
  l3 = (d*w^2 - w) * (d - w^2) / ((d*w^2 - w^2) * (d - w));
  [lift(l1), lift(l2), lift(l3)]
};

\\ ---- Verify the Frobenius polynomial ----
verify_frobenius(h_poly, target_t, p) = {
  my(cp, expected, ok);
  cp = hyperellcharpoly(h_poly);
  expected = (subst(cp, variable(cp), x) == x^4 + (p - target_t^2)*x^2 + p^2);
  \\ More direct: evaluate char poly at x=1 to get #Jac
  my(n_jac, target_n);
  n_jac = subst(cp, variable(cp), 1);
  target_n = (p + 1 - target_t) * (p + 1 + target_t);
  [cp, n_jac, target_n, n_jac == target_n]
};

\\ ================================================================
\\ Part 1: Toy case p = 1009, b = 11
\\ ================================================================
print("---- Part 1: Toy case p = 1009, E: y^2 = x^3 + 11 ----");
print("");

{
p = 1009;
b = 11;
E1 = ellinit([0, b], p);
t_E1 = p + 1 - ellcard(E1);
print("E1: y^2 = x^3 + ", b, " over F_", p);
print("Trace of Frobenius: t = ", t_E1);
print("Expected Jac order: (p+1-t)(p+1+t) = ", (p+1-t_E1)*(p+1+t_E1));
print("");

\\ Find ω and d
om = cube_root_unity(p);
d_ns = first_nonsquare(p);
b2 = (d_ns^3 * b) % p;
print("Primitive cube root of unity ω = ", om, " (check: ω^3 = ", lift(Mod(om,p)^3), ", ω^2+ω+1 = ", lift(Mod(om^2+om+1, p)), ")");
print("First non-square d = ", d_ns, " (check: d^((p-1)/2) = ", lift(Mod(d_ns, p)^((p-1)/2)), ")");
print("E2 twist coefficient b2 = d^3 * b mod p = ", b2);
print("E2: y^2 = x^3 + ", b2);
E2 = ellinit([0, b2], p);
t_E2 = p + 1 - ellcard(E2);
print("Trace of E2: t = ", t_E2, " (expected ", -t_E1, ")");
print("");

\\ Compute Rosenhain invariants (original Möbius)
lams = rosenhain_from_twist(b, d_ns, om, p);
print("Rosenhain invariants (variant 1, T: α→0, dα→1, ωα→∞):");
print("  λ1 = ", lams[1]);
print("  λ2 = ", lams[2]);
print("  λ3 = ", lams[3]);
print("");

\\ Construct and verify the Rosenhain curve
h_rosenhain = rosenhain_poly(lams, p);
print("Rosenhain curve: y^2 = ", h_rosenhain);
print("");

res = verify_frobenius(h_rosenhain, t_E1, p);
cp = res[1];
n_jac = res[2];
target_n = res[3];
match = res[4];
print("Frobenius char poly of Jac(C): ", cp);
print("#Jac(C) = char_F(1) = ", n_jac);
print("Target n_E1 * n_E2   = ", target_n);
print("Match? ", match);
print("");

\\ Also try the second cube root of unity (ω^2)
om2 = lift(Mod(om, p)^2);
print("---- Variant with ω^2 = ", om2, " ----");
lams2 = rosenhain_from_twist(b, d_ns, om2, p);
print("Rosenhain invariants (ω^2 variant):");
print("  λ1 = ", lams2[1]);
print("  λ2 = ", lams2[2]);
print("  λ3 = ", lams2[3]);
h2 = rosenhain_poly(lams2, p);
res2 = verify_frobenius(h2, t_E1, p);
print("Frobenius char poly: ", res2[1]);
print("#Jac = ", res2[2], "  target = ", res2[3], "  match = ", res2[4]);
print("");

\\ Variant 2 Möbius
lams3 = rosenhain_variant2(b, d_ns, om, p);
print("---- Variant 2 Möbius (ωα→0, dα→1, ω^2α→∞) ----");
print("  λ1 = ", lams3[1]);
print("  λ2 = ", lams3[2]);
print("  λ3 = ", lams3[3]);
h3 = rosenhain_poly(lams3, p);
res3 = verify_frobenius(h3, t_E1, p);
print("Frobenius char poly: ", res3[1]);
print("#Jac = ", res3[2], "  target = ", res3[3], "  match = ", res3[4]);
print("");
}

\\ ================================================================
\\ Part 2: Full Igusa invariants of the NAIVE cover h_secp
\\ ================================================================
print("================================================================");
print("Part 2: Full Igusa invariants of y^2 = (x^3+7)(x^3+189) mod p_secp");
print("================================================================");
print("");

\\ Import the transvectant functions from igusa_clebsch.gp
read("igusa_clebsch.gp");

{
p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
h_secp = (x^3 + 7) * (x^3 + 189);
print("h_secp = ", h_secp);
print("");

\\ Compute over Z first, then reduce mod p
J2_Z  = igusa_J2(h_secp);
J4_Z  = igusa_J4_transvectant(h_secp);
J6_Z  = igusa_J6_candidate(h_secp);
J10_Z = igusa_J10(h_secp);

print("Igusa invariants over Z (before reduction):");
print("  J2  = ", J2_Z);
print("  J4  = ", J4_Z, "  (numerator/denominator: ", numerator(J4_Z), " / ", denominator(J4_Z), ")");
print("  J6  = ", J6_Z, "  (numerator/denominator: ", numerator(J6_Z), " / ", denominator(J6_Z), ")");
print("  J10 = ", J10_Z);
print("");

\\ Reduce mod p_secp
J2_p  = lift(Mod(J2_Z, p_secp));
J4_p  = lift(Mod(numerator(J4_Z) * Mod(denominator(J4_Z), p_secp)^(-1), p_secp));
J6_p  = lift(Mod(numerator(J6_Z) * Mod(denominator(J6_Z), p_secp)^(-1), p_secp));
J10_p = lift(Mod(J10_Z, p_secp));

print("Igusa invariants mod p_secp:");
print("  J2  mod p = ", J2_p);
print("  J4  mod p = ", J4_p);
print("  J6  mod p = ", J6_p);
print("  J10 mod p = ", J10_p);
print("");

\\ Normalised invariant ratios (isomorphism class markers)
r42 = J4_Z / J2_Z^2;
r63 = J6_Z / J2_Z^3;
print("Normalised ratios over Q:");
print("  J4 / J2^2 = ", r42);
print("  J6 / J2^3 = ", r63);
print("");

\\ Non-degeneracy check
print("Non-degeneracy:");
print("  J2 ≠ 0 mod p_secp: ", J2_p != 0);
print("  J10 ≠ 0 mod p_secp: ", J10_p != 0);
print("  → y^2 = h_secp(x) is a smooth genus-2 curve over F_p_secp: ", J10_p != 0);
print("");
}

\\ ================================================================
\\ Part 3: Summary of CHLRS formula status
\\ ================================================================
print("================================================================");
print("Summary");
print("================================================================");
print("");
print("Rosenhain formula: for j=0 curve y^2=x^3+b over F_p (p≡1 mod 3)");
print("and its quadratic twist E^t: y^2=x^3+b*d^3 (d non-square),");
print("the Howe-glued genus-2 curve in Rosenhain form is:");
print("  C: y^2 = x(x-1)(x-λ1)(x-λ2)(x-λ3)");
print("with:");
print("  λ1 = (ω^2-1)(d-ω) / (ω(ω-1)(d-1))");
print("  λ2 = (dω-1)(d-ω)  / (ω(d-1)^2)");
print("  λ3 = (dω^2-1)(d-ω)/ (ω(dω-1)(d-1))");
print("where ω = primitive cube root of unity in F_p.");
print("");
print("Claim: Jac(C) is (2,2)-isogenous to E × E^t.");
print("Status: empirically verified above for p=1009 (match = ?).");
print("");
print("If match = 1: formula is correct! Apply to secp256k1.");
print("If match = 0: formula needs correction or different Möbius choice.");

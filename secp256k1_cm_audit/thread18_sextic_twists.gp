\\ thread18_sextic_twists.gp
\\ Check Howe (H1)+(H2)+(H3) for all 15 pairwise combinations
\\ of the 6 sextic twists of secp256k1 (j=0 curve, y^2 = x^3 + 7).
\\
\\ Run: gp --stacksize 128000000 -q thread18_sextic_twists.gp

default(parisize, 128000000);
default(timer, 0);

\\ ====================================================================
\\ Utilities
\\ ====================================================================

\\ Return degree-pattern vector from factored polynomial over F_p
deg_pattern(f) = {
  my(nf = matsize(f)[1], v = []);
  for (r = 1, nf, v = concat(v, poldegree(f[r,1])));
  Vec(v)
};

\\ Howe H2: same 2-torsion Galois module structure (same degree pattern for x^3+b)
same_2torsion(b1, b2, p) = {
  my(f1 = factor(Mod(x^3 + b1, p)));
  my(f2 = factor(Mod(x^3 + b2, p)));
  deg_pattern(f1) == deg_pattern(f2)
};

\\ Check all 3 Howe conditions for a pair (ni, nj, bi, bj) over prime p
\\ ni, nj = group orders; bi, bj = curve constants (y^2 = x^3 + b)
howe_check(ni, nj, bi, bj, p) = {
  my(H1, H2, H3, g12, all3);
  H1 = (ni != nj);
  H2 = same_2torsion(bi, bj, p);
  g12 = gcd(ni, nj);
  H3 = (g12 == 1);
  all3 = H1 && H2 && H3;
  [H1, H2, H3, g12, all3]
};

\\ ====================================================================
\\ Part A: Toy prime (p0 = 43, p0 mod 6 = 1), explicit group orders
\\ ====================================================================
print("================================================================");
print("Thread 18: Howe gluing test for secp256k1 sextic twists");
print("================================================================");
print("");
print("=== Part A: toy prime p0 = 43 ===");
print("");

p0 = 43;
printf("p0 = %d, p0 %% 6 = %d\n", p0, p0 % 6);

\\ primitive root mod p0
g0 = lift(znprimroot(p0));
printf("primitive root g0 = %d\n\n", g0);

\\ Group orders of the 6 sextic twists E_k: y^2 = x^3 + g0^k
orders0 = vector(6);
bs0 = vector(6);
{
  for (k = 0, 5,
    bk = lift(Mod(g0, p0)^k);
    bs0[k+1] = bk;
    Ek = ellinit([0, 0, 0, 0, bk], p0);
    orders0[k+1] = ellcard(Ek);
    printf("  k=%d: b=%d, #E_%d = %d, trace = %d\n",
      k, bk, k, orders0[k+1], p0+1-orders0[k+1]);
  );
}
print("");

\\ 2-torsion splitting
print("2-torsion degree patterns (x^3 + g0^k mod p0):");
{
  for (k = 0, 5,
    bk = bs0[k+1];
    f = factor(Mod(x^3 + bk, p0));
    dp = deg_pattern(f);
    printf("  k=%d: b=%d => %s\n", k, bk, Vec(dp));
  );
}
print("");

\\ All 15 pairs
print("All 15 pairs (p0=43):");
printf("  %-10s  H1  H2  H3  gcd  result\n", "(E_i, E_j)");
hits0 = 0;
{
  for (i = 0, 4,
    for (j = i+1, 5,
      r = howe_check(orders0[i+1], orders0[j+1], bs0[i+1], bs0[j+1], p0);
      if (r[5], hits0++);
      printf("  (E_%d, E_%d)     %d   %d   %d   %d    %s\n",
        i, j, r[1], r[2], r[3], r[4], if(r[5], "GLUE", "--"));
    );
  );
}
printf("\nPairs with H1+H2+H3 at p0=%d: %d / 15\n\n", p0, hits0);

\\ ====================================================================
\\ Part B: secp256k1 prime — group orders via CM (Eisenstein) theory
\\ ====================================================================
print("=== Part B: secp256k1 prime (CM theory) ===");
print("");

p256k1 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n256k1 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t256k1 = p256k1 + 1 - n256k1;

printf("p = (256-bit prime)\n");
printf("n = #E_0(F_p)  [secp256k1, prime order]\n");
printf("t = p+1-n = %d\n", t256k1);
printf("p mod 6 = %d  (6 | p-1, so exactly 6 sextic twists exist)\n\n", p256k1 % 6);

\\ secp256k1 CM: Z[omega], omega = e^{2pi i/3}, norm form a^2 - ab + b^2 = p.
\\ Frobenius pi = a + b*omega; trace = 2a - b = t.
\\ b = 2a - t; substitute into norm:
\\   a^2 - a(2a-t) + (2a-t)^2 = p
\\   3a^2 - 3at + t^2 - p = 0
\\   discriminant: 12p - 3t^2
Delta = 12*p256k1 - 3*t256k1^2;
sq = sqrtint(abs(Delta));
is_sq = (sq^2 == Delta);
printf("12p - 3t^2 = %d  (>0: %d, perfect square: %d)\n", Delta, Delta > 0, is_sq);

if (!is_sq,
  error("CM parametrization: 12p-3t^2 is not a perfect square");
);

\\ Solve for (a, b)
a_cm = (3*t256k1 + sq) \ 6;
b_cm = 2*a_cm - t256k1;
ok1 = (a_cm^2 - a_cm*b_cm + b_cm^2 == p256k1) && (2*a_cm - b_cm == t256k1);
if (!ok1,
  a_cm = (3*t256k1 - sq) \ 6;
  b_cm = 2*a_cm - t256k1;
);
printf("CM parametrization: a = %d, b = %d\n", a_cm, b_cm);
printf("  Verify: a^2 - ab + b^2 == p? %d\n", a_cm^2 - a_cm*b_cm + b_cm^2 == p256k1);
printf("  Verify: 2a - b == t? %d\n\n", 2*a_cm - b_cm == t256k1);

\\ The 6 sextic twist traces are Re(2 * pi * zeta6^k) where zeta6 = unit of Z[omega].
\\ Units of Z[omega]: {1, omega, omega^2, -1, -omega, -omega^2}.
\\ pi = a + b*omega; trace(pi*u_k) = pi*u_k + conjugate(pi*u_k).
\\
\\ Explicit formula (from Eisenstein arithmetic):
\\   k=0 (u=1):      trace = 2a - b
\\   k=1 (u=omega):  trace = -a - b
\\   k=2 (u=omega^2):trace = 2b - a
\\   k=3 (u=-1):     trace = b - 2a       [= quadratic twist of k=0]
\\   k=4 (u=-omega): trace = a + b         [= quadratic twist of k=1]
\\   k=5 (u=-omega^2):trace = a - 2b      [= quadratic twist of k=2]

a = a_cm; b = b_cm;
traces = [2*a - b, -a - b, 2*b - a, b - 2*a, a + b, a - 2*b];
orders = vector(6);
print("6 sextic twist group orders (from CM theory):");
{
  for (k = 0, 5,
    orders[k+1] = p256k1 + 1 - traces[k+1];
    printf("  k=%d: trace=%d, #E_%d = p+1-trace\n", k, traces[k+1], k);
    \\ Factor the order for display
    printf("       #E_%d mod 4 = %d, #E_%d mod 2 = %d\n",
      k, orders[k+1] % 4, k, orders[k+1] % 2);
  );
}
printf("\nSanity: k=0 trace == t? %d\n", traces[1] == t256k1);
printf("Sanity: k=0 order == n? %d\n\n", orders[1] == n256k1);

\\ 2-torsion type from group order:
\\ - If 4 | #E_k: x^3 + 7*g^k splits (all 3 non-trivial 2-torsion F_p-rational)
\\ - If 2 | #E_k but 4 ∤: partial split (1 F_p-rational 2-torsion point)
\\ - If #E_k odd: x^3 + 7*g^k irreducible (no F_p-rational non-trivial 2-torsion)
\\
\\ For j=0 with p ≡ 1 mod 3: the 2-torsion poly x^3+7*g^k is either
\\   irreducible (when -7*g^k is not a cube mod p) or splits into 3 linear factors.
\\   There is no "1 linear + 1 quadratic" case for j=0 over F_p with p≡1(mod3).
\\ Criterion: "splits" iff 4 | #E_k (all non-trivial 2-torsion is F_p-rational).

print("2-torsion type of each twist:");
split_ks = []; irred_ks = [];
{
  for (k = 0, 5,
    ok = orders[k+1] % 4;
    if (ok == 0,
      ttype = "split (4 | #E)";
      split_ks = concat(split_ks, [k])
    ,
      ttype = Str("irred (#E mod 4 = ", ok, ")");
      irred_ks = concat(irred_ks, [k])
    );
    printf("  k=%d: type = %s\n", k, ttype);
  );
}
printf("\n'Split' k's: %s\n", Vec(split_ks));
printf("'Irred' k's: %s\n\n", Vec(irred_ks));

\\ H2 classification:
\\ - Pairs within split group: H2 holds
\\ - Pairs within irred group: H2 holds
\\ - Cross-group pairs: H2 FAILS (different Galois module structure)
\\ Expected: 2 split + 4 irred => 1 + 6 = 7 same-type pairs; 8 cross pairs.
\\
\\ For H2, we use the group-order criterion (proxy for actual poly factoring).
\\ This is valid for j=0 / p≡1(mod 3): no intermediate "partial split" case.

print("All 15 pairs (secp256k1, from CM theory):");
printf("  %-10s  H1  H2  H3  gcd  result\n", "(E_i,E_j)");
hits = 0;
fails_H1 = []; fails_H2 = []; fails_H3 = [];
{
  for (i = 0, 4,
    for (j = i+1, 5,
      ni = orders[i+1]; nj = orders[j+1];
      H1 = (ni != nj);
      \\ H2 via group-order proxy: same_type iff both 4|order or both 4∤order
      H2 = ((ni % 4 == 0) == (nj % 4 == 0));
      g12 = gcd(ni, nj);
      H3 = (g12 == 1);
      all3 = H1 && H2 && H3;
      if (all3, hits++);
      if (!H1, fails_H1 = concat(fails_H1, [[i,j]]));
      if (!H2, fails_H2 = concat(fails_H2, [[i,j]]));
      if (!H3, fails_H3 = concat(fails_H3, [[i,j]]));
      \\ Show gcd only if small for readability
      gcd_str = if(g12 == 1, "1", Str(g12));
      printf("  (E_%d,E_%d)     %d   %d   %d   %s    %s\n",
        i, j, H1, H2, H3, gcd_str, if(all3, "GLUE", "--"));
    );
  );
}
printf("\nPairs passing H1+H2+H3 (secp256k1): %d / 15\n", hits);
printf("H1 failures: %d\n", #fails_H1);
printf("H2 failures: %d\n", #fails_H2);
printf("H3 failures: %d\n\n", #fails_H3);

\\ ====================================================================
\\ Summary
\\ ====================================================================
print("================================================================");
print("Summary");
print("================================================================");
print("");
print("secp256k1 has 6 sextic twists E_0..E_5 (j=0, p ≡ 1 mod 6).");
print("2-torsion splits: 2 twists have x^3+7*g^k splitting ('split'),");
print("                  4 twists have x^3+7*g^k irreducible ('irred').");
print("H2 holds for:   1 pair (split-split) + 6 pairs (irred-irred) = 7 pairs.");
print("H2 fails for:   8 cross-group pairs.");
print("");
print("Among the 7 H2-passing pairs:");
print("  H1 (non-isogenous) — fails only if two distinct twists have the");
print("    same group order (quadratic twists of each other with t=0 would");
print("    fail, but t ≠ 0 for secp256k1, so H1 holds for all 7).");
print("  H3 (gcd=1) — E_0 has prime order n, guaranteeing gcd=1 for all");
print("    pairs involving E_0; for other pairs must check numerics.");
print("");
printf("RESULT: %d of 15 pairs admit (2,2)-Howe gluing to a genus-2 Jacobian.\n", hits);
print("");
print("ECDLP implication:");
print("  Each gluing pair yields a genus-2 curve C/F_p with Jac(C) -> E_i x E_j.");
print("  Best HCDLP on genus-2 over F_p: Pollard rho, cost O(p).");
print("  ECDLP cost: O(sqrt(n)) ≈ O(p^{1/2}).  No cryptographic speedup.");
print("");
print("STRUCTURAL OBSERVATION:");
print("  The 'irred' quadratic-twist pairs (E_k, E_{k+3}) with k in irred-group");
print("  each have #E_k * #E_{k+3} ≈ p^2 but with gcd(#E_k, #E_{k+3}) > 1");
print("  in general (both divisible by the same small factors from p+1±t).");
print("  These are the H3-failure candidates; see gcd column above.");
print("");
print("Thread 18 complete.");

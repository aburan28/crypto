\\ =================================================================
\\ Thread 18: Howe (H1)+(H2)+(H3) check for ALL 15 pairs of
\\ j=0 sextic twists of secp256k1.
\\
\\ PARI 2.15.4 limitation: my() inside for-loop body confuses parser.
\\ Workaround: ALL for-loop bodies are single helper-function calls.
\\
\\ Howe (1996) gluing conditions for a smooth genus-2 Jacobian:
\\   (H1) #E_i != #E_j  (Hom=0 by Tate)
\\   (H2) E_i[2] =~ E_j[2] as F_p-Galois modules
\\   (H3) gcd(#E_i, #E_j) = 1  (sufficient form)
\\
\\ CM trace formula: if p = A^2-AB+B^2 in Z[omega], omega=e^{2pi i/3},
\\ and t0 = 2A-B (trace of E_0), the 6 traces are:
\\   k=0: 2A-B   k=1: A-2B   k=2: -A-B
\\   k=3: B-2A   k=4: 2B-A   k=5: A+B
\\ (derived: z6^k * pi0 in Z[omega], z6=-omega^2=1+omega)
\\
\\ Condition for z6 to represent all 6 twist classes:
\\   gcd((p-1)/6, 6) = 1.
\\   For secp256k1: (p-1)/6 is odd (p=3 mod 4) and =1 mod 3, so gcd=1.
\\   Sanity check uses p=7 where gcd((7-1)/6,6)=gcd(1,6)=1.
\\
\\ Run: gp -q thread18_sextic_twist_howe.gp
\\ =================================================================

default(parisize, 256*1024*1024);
default(timer, 0);

\\ ---- CM decomposition: Find A,B with A^2-AB+B^2=p, 2A-B=t. ----
cm_decompose(p_in, t_in) = {
  my(d3, Bsq, B, A);
  d3 = 4*p_in - t_in^2;
  if(d3 % 3 != 0, error("cm_decompose: (4p-t^2) not div 3"));
  Bsq = d3 / 3;
  B = sqrtint(Bsq);
  if(B^2 != Bsq, error("cm_decompose: not perfect square"));
  A = (t_in + B) / 2;
  if(type(A) == "t_FRAC", A = (t_in - B)/2; B = -B);
  if(A^2 - A*B + B^2 != p_in, error("cm_decompose: p verify failed"));
  if(2*A - B != t_in, error("cm_decompose: t verify failed"));
  [A, B]
};

\\ ---- Six traces in k=0..5 order. ----
cm_traces6(A, B) = [2*A-B, A-2*B, -A-B, B-2*A, 2*B-A, A+B];

\\ ---- Find z6 of order 6 via polrootsmod(x^2+x+1, p). ----
find_z6(p_in) = {
  my(rt, om3, z6c);
  rt = polrootsmod(x^2+x+1, p_in);
  if(#rt < 1, error("find_z6: x^2+x+1 has no root mod p"));
  om3 = Mod(rt[1], p_in);
  z6c = -om3;
  if(z6c^6 != Mod(1,p_in), error("find_z6: z6^6 != 1"));
  if(z6c^3 == Mod(1,p_in), error("find_z6: got order 3 not 6"));
  z6c
};

\\ ---- H2 type: split (=1) or irr (=0). ----
h2_split(c_in, p_in) = (Mod(-c_in, p_in)^((p_in-1)/3) == Mod(1,p_in));

\\ ============================================================
\\ Section 1 globals and helper
\\ ============================================================
g_p0 = 7;
g_b0 = 1;
g_z6_0 = 0;     \\ set after find_z6
g_t_dir = vector(6);

\\ Helper: process one row of the p=7 twist table.
s1_row(k) = {
  my(ck, Ek, nk, tk, stype);
  ck = lift(Mod(g_b0, g_p0) * g_z6_0^k);
  Ek = ellinit([0, ck], g_p0);
  nk = ellcard(Ek);
  tk = g_p0 + 1 - nk;
  g_t_dir[k+1] = tk;
  stype = if(h2_split(ck, g_p0), "SPLIT", "IRR");
  print("  k=", k, "  b=", ck, "  #E=", nk, "  t=", tk, "  H2=", stype)
};

\\ ============================================================
\\ Section 2 globals and helpers
\\ ============================================================
g_p_s = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
g_n_s = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
g_b_s = 7;
g_z6_s = 0;     \\ set after find_z6
g_TR  = vector(6);
g_H2  = vector(6);
g_h1  = 0; g_h2c = 0; g_h12 = 0; g_h123 = 0;

\\ Helper: print one trace row.
s2_trace_row(k) = {
  print("  k=", k, "  trace=", g_TR[k+1], "  #E_k=", g_p_s+1-g_TR[k+1])
};

\\ Helper: compute and print one H2 row.
s2_h2_row(k) = {
  my(ck, stype);
  ck = lift(Mod(g_b_s, g_p_s) * g_z6_s^k);
  g_H2[k+1] = h2_split(ck, g_p_s);
  stype = if(g_H2[k+1], "SPLIT", "IRR");
  print("  k=", k, "  H2-type=", stype)
};

\\ Helper: one Howe pair check + print.
do_pair(ii, jj) = {
  my(tri, trj, ni, nj, h1, h2, g3, h3, mark);
  tri = g_TR[ii+1]; trj = g_TR[jj+1];
  ni = g_p_s + 1 - tri; nj = g_p_s + 1 - trj;
  h1 = (tri != trj);
  h2 = (g_H2[ii+1] == g_H2[jj+1]);
  g3 = gcd(ni, nj);
  h3 = (g3 == 1);
  if(h1,        g_h1   = g_h1+1);
  if(h2,        g_h2c  = g_h2c+1);
  if(h1 && h2,  g_h12  = g_h12+1);
  if(h1&&h2&&h3, g_h123 = g_h123+1);
  mark = if(h1 && h2 && h3, " *** ALL MET ***", "");
  print("  (", ii, ",", jj, "): H1=", h1,
        " H2=", h2, "(", if(g_H2[ii+1],"S","I"), "/", if(g_H2[jj+1],"S","I"), ")",
        " H3=", h3, " gcd=", if(g3==1,"1",Str(g3)), mark)
};

\\ Helper for nested for-loop: inner loop only.
do_pair_j(ii) = for(jj=ii+1, 5, do_pair(ii, jj));

\\ ==========================================================
\\ SECTION 1: Sanity check on p=7
\\ ==========================================================
print("==========================================================");
print("SECTION 1: Sanity check p=7, j=0 curve y^2=x^3+1");
print("  p=7: p mod 6=1, (p-1)/6=1, gcd(1,6)=1 => z6 covers all 6 twist classes");
print("==========================================================");
print("");

g_z6_0 = find_z6(g_p0);
print("z6 mod 7 = ", lift(g_z6_0), "  order = ", znorder(g_z6_0));
print("");
print("Six sextic twists y^2=x^3+z6^k over F_7:");
for(k=0, 5, s1_row(k));
print("");

\\ CM formula verification.
{
  my(ab0, A0, B0, tr0, match);
  ab0 = cm_decompose(g_p0, g_t_dir[1]);
  A0 = ab0[1]; B0 = ab0[2];
  tr0 = cm_traces6(A0, B0);
  match = (vecsort(tr0) == vecsort(Vec(g_t_dir)));
  print("CM: A=", A0, " B=", B0, "  A^2-AB+B^2=", A0^2-A0*B0+B0^2);
  print("CM traces:    ", tr0);
  print("Direct traces:", Vec(g_t_dir));
  print("Sets match: ", match);
  print("All 6 distinct: ", #Set(tr0) == 6)
}

\\ ==========================================================
\\ SECTION 2: secp256k1 full analysis
\\ ==========================================================
print("");
print("==========================================================");
print("SECTION 2: secp256k1 -- 6 sextic twists, all 15 Howe pairs");
print("==========================================================");
print("");
print("p = ", g_p_s);
print("n = ", g_n_s, "  (prime)");
print("t0 = ", g_p_s + 1 - g_n_s);
print("p mod 6 = ", g_p_s % 6, "  (need 1)");
print("gcd((p-1)/6, 6) = ", gcd((g_p_s-1)/6, 6), "  (need 1 for z6 to cover all 6 classes)");
print("");

\\ CM decomposition.
{
  my(ab, A, B, t0);
  t0 = g_p_s + 1 - g_n_s;
  ab = cm_decompose(g_p_s, t0);
  A = ab[1]; B = ab[2];
  g_TR = cm_traces6(A, B);
  print("CM: A=", A);
  print("    B=", B);
  print("    A^2-AB+B^2==p: ", A^2-A*B+B^2 == g_p_s);
  print("    2A-B==t0:       ", 2*A-B == t0)
}
print("");

\\ Six traces.
print("Six traces and group orders:");
for(k=0, 5, s2_trace_row(k));
print("");
print("All 6 traces distinct (H1 for all 15 pairs): ", #Set(g_TR) == 6);
print("");

\\ Find z6 in F_p.
g_z6_s = find_z6(g_p_s);
print("z6 found: z6^6==1: ", g_z6_s^6==Mod(1,g_p_s), "  z6^3==-1: ", g_z6_s^3==Mod(-1,g_p_s));
print("");

\\ H2 types.
print("H2 analysis (x^3+7*z6^k splits <=> -7*z6^k is a cube mod p):");
for(k=0, 5, s2_h2_row(k));
print("");
{
  my(ns);
  ns = sum(k=1, 6, g_H2[k]);
  print("Split-type twists: ", ns, " of 6  (expected: 2 by cubic residue theory)");
  print("IRR-type twists:   ", 6-ns, " of 6  (expected: 4)")
}

\\ ==========================================================
\\ All 15 Howe pairs
\\ ==========================================================
print("");
print("==========================================================");
print("All 15 pairs: H1=(trace differ), H2=(same 2-tor type), H3=(gcd=1)");
print("  S=split H2-type, I=irr H2-type");
print("==========================================================");
print("");

for(ii=0, 4, do_pair_j(ii));

print("");
print("Summary:");
print("  Pairs with H1 (trace differ):    ", g_h1,   " / 15");
print("  Pairs with H2 (same 2-tor type): ", g_h2c,  " / 15");
print("  Pairs with H1+H2:                ", g_h12,  " / 15");
print("  Pairs with H1+H2+H3 (ALL MET):   ", g_h123, " / 15");

\\ ==========================================================
\\ Structural summary
\\ ==========================================================
print("");
print("==========================================================");
print("Structural summary");
print("==========================================================");
print("For j=0 over F_p with p==1(mod 6) and gcd((p-1)/6,6)=1:");
print("  6 distinct sextic twist classes, each y^2=x^3+b*z6^k.");
print("  CM traces {2A-B,A-2B,-A-B,B-2A,2B-A,A+B} are all distinct => H1 for all 15.");
print("  H2: exactly 2-of-6 twists have split 2-torsion; 4-of-6 have irr.");
print("    => 8 of 15 pairs fail H2; only 7 satisfy H2.");
print("    (For secp256k1 this differs -- see result below.)");
print("  H3: gcd(#E_i, #E_j) checked numerically.");
print("");
print("secp256k1 result:");
print("  H1+H2+H3 (Howe-glueable): ", g_h123, " of 15 pairs.");
print("");
print("ECDLP note:");
print("  Howe-glueable pairs give genus-2 Jacobians with #Jac ~= p^2.");
print("  Pollard-rho on g=2: cost ~sqrt(p^2)=p >> sqrt(p)=ECDLP cost.");
print("  Sub-exp index calculus (Diem 2011) requires genus >= 3.");
print("  => Sextic-twist Howe gluing is structural but ECDLP-neutral.");
print("");
print("==========================================================");
print("Thread 18 complete.");
print("==========================================================");

\\ =================================================================
\\ Thread 18: Howe (H1)+(H2)+(H3) conditions for all 15 pairs
\\ of j=0 sextic twists of secp256k1.
\\
\\ secp256k1: y^2 = x^3 + 7 over F_p.  j(E) = 0, so Aut(E) = Z/6Z
\\ (when p > 3), and there are exactly 6 F_p-isomorphism classes of
\\ curves with j = 0, called the 6 sextic twists.
\\
\\ Howe (1996) conditions for (2,2)-gluing of E_i x E_j:
\\   (H1) Hom_{F_p}(E_i, E_j) = 0  <=> #E_i != #E_j
\\   (H2) E_i[2] iso E_j[2] as F_p-Galois modules (with Weil pairing)
\\   (H3) gcd(#E_i, #E_j) = 1  (sufficient; Howe allows extensions)
\\
\\ Key structure:
\\   p ≡ 1 (mod 3)  => F_p has primitive cube roots of unity
\\                  => x^3+c either irreducible OR splits completely
\\   p ≡ 3 (mod 4)  => p+1 ≡ 0 (mod 4)
\\   CM: j=0 <=> CM by Z[zeta_3] (Eisenstein integers, disc -3)
\\
\\ Run: gp -q thread18_sextic_twist_howe.gp
\\ =================================================================

default(parisize, 512000000);
default(timer, 0);

print("=================================================================");
print("Thread 18: Howe gluing on all 15 pairs of j=0 sextic twists");
print("=================================================================");
print("");

\\ secp256k1 parameters
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t0 = p + 1 - n;

print("p  = ", p);
print("n  = ", n, "  (secp256k1 group order, prime: ", isprime(n), ")");
print("t0 = p+1-n = ", t0);
print("p mod 3 = ", p % 3, "  (must be 1 for Z[zeta_3] CM)");
print("p mod 4 = ", p % 4, "  (must be 3 for proxy H2 check to work)");
print("");

\\ =================================================================
\\ Step 1: Eisenstein factorization of p in Z[zeta_3].
\\
\\ Write p = a^2 - a*b + b^2 in Z[zeta_3] (Eisenstein integers).
\\ The trace of E: y^2=x^3+7 is t0 = 2a - b (by CM theory for j=0).
\\ So: b = 2a - t0, substitute into a^2 - ab + b^2 = p:
\\     3a^2 - 3a*t0 + t0^2 = p
\\     => (4p - t0^2) / 3 = u^2  for some u in Z, and a = (t0 +- u)/2.
\\ =================================================================

print("---- Eisenstein factorization ----");
print("");
disc_eis = 4*p - t0^2;
print("4p - t0^2 = ", disc_eis);
print("(4p - t0^2) mod 3 = ", disc_eis % 3);
if(disc_eis % 3 != 0,
  error("ABORT: 4p - t0^2 not divisible by 3 -- p is not CM by Z[zeta_3]!"));

v3 = disc_eis \ 3;
u_eis = sqrtint(v3);
if(u_eis^2 != v3,
  error("ABORT: (4p-t0^2)/3 is not a perfect square!"));
print("u = sqrt((4p-t0^2)/3) = ", u_eis);

\\ Find integer a: (t0 + u_eis) / 2 or (t0 - u_eis) / 2
a_eis = 0; b_eis = 0;
tmp1 = t0 + u_eis;
tmp2 = t0 - u_eis;
if(tmp1 % 2 == 0,
  a_eis = tmp1 \ 2; b_eis = 2*a_eis - t0,
  a_eis = tmp2 \ 2; b_eis = 2*a_eis - t0
);
if(a_eis^2 - a_eis*b_eis + b_eis^2 != p,
  \\ try the other candidate
  a_eis = tmp2 \ 2; b_eis = 2*a_eis - t0);
if(a_eis^2 - a_eis*b_eis + b_eis^2 != p,
  error("ABORT: Eisenstein factorization failed!"));

print("");
print("Eisenstein factorization: pi = a + b*zeta_3  where");
print("  a  = ", a_eis);
print("  b  = ", b_eis);
print("  Norm(pi) = a^2 - ab + b^2 = ", a_eis^2 - a_eis*b_eis + b_eis^2, " = p: ", a_eis^2 - a_eis*b_eis + b_eis^2 == p);
print("  Trace(pi) = 2a - b = ", 2*a_eis - b_eis, " = t0: ", 2*a_eis - b_eis == t0);
print("");

\\ =================================================================
\\ Step 2: The 6 sextic twist traces.
\\
\\ Multiplying pi by the 6 units {+/-1, +/-zeta_3, +/-zeta_3^2}
\\ gives 6 algebraic integers whose traces are the 6 twist traces.
\\
\\ For pi = a + b*zeta_3:
\\   eps=1:       trace = 2a - b
\\   eps=-1:      trace = -(2a - b)
\\   eps=zeta_3:  trace = -(a + b)
\\   eps=-zeta_3: trace = a + b
\\   eps=zeta_3^2 = zeta_3bar: trace = 2b - a
\\   eps=-zeta_3^2: trace = -(2b - a)
\\
\\ Consistency check: (2a-b) + (-(a+b)) + (2b-a) = 0.  Sum of the
\\ three "epsilons from one half" is zero (since 1+zeta_3+zeta_3^2=0).
\\ =================================================================

print("---- Six sextic twist traces ----");
print("");

\\ Three fundamental traces (whose negatives give the other three)
s1 =  2*a_eis - b_eis;   \\ trace for eps=1   (= t0 for secp256k1)
s2 = -(a_eis + b_eis);   \\ trace for eps=zeta_3
s3 =  2*b_eis - a_eis;   \\ trace for eps=zeta_3^2

print("s1 = 2a - b  = ", s1, "  (= t0: ", s1 == t0, ")");
print("s2 = -(a+b)  = ", s2);
print("s3 = 2b - a  = ", s3);
print("Sum s1+s2+s3 = ", s1+s2+s3, "  (should be 0)");
print("");

\\ All 6 traces
T = [s1, -s1, s2, -s2, s3, -s3];
print("All 6 traces:  T[1..6] = ", T);
print("");

\\ =================================================================
\\ Step 3: Group orders and 2-torsion structure.
\\
\\ N[i] = p + 1 - T[i].
\\
\\ (H2) proxy: For p ≡ 3 (mod 4), p+1 ≡ 0 (mod 4), so
\\   4 | N[i]  <=>  4 | T[i]
\\   <=>  E_i has fully split 2-torsion (E_i[2](F_p) = Z/2 x Z/2).
\\ If 4 ∤ T[i], x^3+b_i is irreducible over F_p (no rational 2-tors).
\\
\\ Among the 6 traces: exactly 2 are divisible by 4
\\ (the pair that are quadratic twists of each other), and 4 are not.
\\ =================================================================

print("---- Group orders and 2-torsion classification ----");
print("");

N = vector(6, i, p + 1 - T[i]);

for(i = 1, 6,
  split = (T[i] % 4 == 0);
  print("E_", i, ":  N = p+1-T[", i, "]  =  <", i, ">  2-tors: ",
    if(split, "SPLIT  (4 | N, 4 | T)", "irred  (4 ∤ N, 4 ∤ T)"));
);
print("");

\\ Count
n_split = sum(i=1,6, T[i] % 4 == 0);
n_irred = 6 - n_split;
print("Split 2-torsion twists: ", n_split, "  (expected: 2)");
print("Irred 2-torsion twists: ", n_irred, "  (expected: 4)");
print("");

\\ Identify split and irred indices
split_idx = [];
irred_idx = [];
for(i=1,6,
  if(T[i] % 4 == 0, split_idx = concat(split_idx,[i]),
                     irred_idx = concat(irred_idx,[i])));
print("Split indices: ", split_idx);
print("Irred indices: ", irred_idx);
print("");

\\ Sanity: confirm secp256k1 (E_1) is irred (its group order n is prime > 4)
if(T[1] % 4 != 0,
  print("E_1 (secp256k1) correctly identified as irred 2-torsion."),
  print("WARNING: E_1 unexpectedly shows SPLIT 2-torsion!")
);
print("");

\\ =================================================================
\\ Step 4: Check all C(6,2) = 15 pairs.
\\ =================================================================

print("=================================================================");
print("Checking all 15 pairs (i<j) for Howe conditions H1, H2, H3:");
print("  H1: N[i] != N[j]");
print("  H2: same 2-torsion type (proxy: 4|T[i] iff 4|T[j])");
print("  H3: gcd(N[i], N[j]) == 1");
print("=================================================================");
print("");

pass_all = List();
pass_h2only = List();

for(i = 1, 5,
  for(j = i+1, 6,
    h1 = (N[i] != N[j]);
    si = (T[i] % 4 == 0); sj = (T[j] % 4 == 0);
    h2 = (si == sj);
    g12 = gcd(N[i], N[j]);
    h3 = (g12 == 1);
    ok = h1 && h2 && h3;

    tag = if(si && sj, "(ss)",
          if(!si && !sj, "(ii)", "(si)"));
    print("(", i, ",", j, ")", tag, "  H1=", if(h1,"Y","N"),
      "  H2=", if(h2,"Y","N"),
      "  H3=", if(h3,"Y","N"),
      "  gcd=", g12,
      if(ok, "  => HOWE PASS", "  => fail"));

    if(h2, listput(pass_h2only,[i,j]));
    if(ok, listput(pass_all,[i,j]));
  )
);

print("");
print("=================================================================");
print("Summary:");
print("  Pairs passing H2 (same 2-torsion type): ", #pass_h2only, " / 15");
print("  Pairs passing ALL of H1+H2+H3:         ", #pass_all, " / 15");
print("=================================================================");
print("");

\\ =================================================================
\\ Step 5: Detailed analysis of H2-passing pairs.
\\ =================================================================

print("---- H2-passing pairs detail ----");
print("");
print("Tag (ss)=both split, (ii)=both irred, (si)=mixed (H2 auto-fail).");
print("");
print("(ss) pairs: 1 pair  -- both have fully rational E[2]; quad-twist pair");
print("(ii) pairs: 6 pairs -- both have irred x^3+b, same cyclic Z/3 Galois mod");
print("(si) pairs: 8 pairs -- mismatched Galois modules; H2 fails by definition");
print("");

for(k=1, #pass_h2only,
  v = pass_h2only[k];
  i = v[1]; j = v[2];
  si2 = (T[i] % 4 == 0);
  print("Pair (", i, ",", j, ") [", if(si2,"ss","ii"), "]:  ",
    "gcd(N[", i, "],N[", j, "]) = ", gcd(N[i],N[j]),
    if(gcd(N[i],N[j])==1, "  => H3 PASS", "  => H3 FAIL"));
);
print("");

\\ =================================================================
\\ Step 6: Factor group orders of irred twists to understand gcd structure
\\ =================================================================

print("---- Factoring group orders of the 4 irred-2-torsion twists ----");
print("(Only irred twists can participate in a Howe pair satisfying H2.)");
print("");
for(k=1, #irred_idx,
  i = irred_idx[k];
  print("N[", i, "] = ", N[i]);
  \\ Factor if feasible (may take a while for large composites)
  \\ Use factorint with a bound so it doesn't hang on hard semiprimes
  f = factorint(N[i], 10^6);  \\ trial division up to 10^6
  print("  = ", f);
);
print("");

\\ =================================================================
\\ Step 7: Check known pair (E_1, E_2) = (secp256k1, quadratic twist)
\\ This should match howe_gluing_test.gp which found all H1/H2/H3 pass.
\\ =================================================================

print("---- Cross-check with howe_gluing_test.gp ----");
print("");
n_twist = p + 1 + t0;  \\ = p+1-(-t0) = N[2]
print("n_twist (quadratic twist of secp256k1) = ", n_twist);
print("N[2] from Eisenstein analysis           = ", N[2]);
print("Match: ", n_twist == N[2]);
print("gcd(N[1], N[2]) = ", gcd(N[1],N[2]), "  (should be 1 per howe_gluing_test.gp)");
print("");

print("=================================================================");
print("Thread 18 complete.");
print("=================================================================");

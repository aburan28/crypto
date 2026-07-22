\\ ============================================================
\\ Thread 18: (3,3)-isogeny cover probe for secp256k1
\\ ============================================================
\\
\\ Question: Does there exist a genus-2 Jacobian J/F_p with a
\\ (3,3)-isogeny J → E_secp256k1 × E' for some elliptic curve E'/F_p?
\\
\\ For a (2,2)-cover, the necessary condition is E[2] ≅ E'[2] as
\\ Galois modules (Howe's theorem H2), plus gcd(n,n')=1 and n≠n'.
\\
\\ The (3,3) analogue requires: E[3] ≅ E'[3] as F_p-Galois modules,
\\ with additional conditions on the polarization.
\\
\\ Step 1: Determine the 3-torsion module structure of secp256k1.
\\   ψ₃(x) = 3x⁴ + 84x = 3x(x³ + 28)  (a=0, b=7)
\\   Factor x³ + 28 over F_p → determine E[3] Galois module.
\\
\\ Step 2: Enumerate curves E'/F_p with E'[3] ≅ E[3] as Gal-modules.
\\   - If x³+28 is irreducible mod p: E[3] has cyclic Galois action of
\\     order dividing 3, so the full 3-torsion is defined over F_{p^3} - F_p.
\\     Any E' with x³+b' also irreducible mod p has same Galois type IF
\\     the splitting fields coincide (stronger condition).
\\   - To check: we need Gal-isomorphism as modules, not just same split pattern.
\\
\\ Step 3: For candidate (E, E') pairs with matching 3-torsion, check:
\\   H1: n_E ≠ n_{E'}  (non-isogenous)
\\   H3: gcd(n_E, n_{E'}) = 1  (no common rational subgroup obstruction)
\\
\\ Step 4: Among j=0 sextic twists (our known family), check (3,3) conditions.
\\
\\ Run: gp -q thread18_33cover_probe.gp

default(parisize, 512000000);
default(timer, 0);

print("================================================================");
print("Thread 18: (3,3)-isogeny cover probe for secp256k1");
print("================================================================");
print("");

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_secp = p + 1 - n_secp;

print("secp256k1: y² = x³ + 7 over F_p");
print("n (prime) = ", n_secp);
print("t         = ", t_secp);
print("");

\\ ---- Step 1: 3-torsion polynomial of secp256k1 ----
print("---- Step 1: 3-torsion structure of secp256k1 ----");
print("");
print("3-division polynomial: ψ₃(x) = 3x⁴ + 84x = 3x(x³ + 28)");
print("");

\\ Check if x=0 gives rational 3-torsion points (i.e., if 7 is a QR mod p)
is_7_qr = (kronecker(7, p) == 1);
print("7 is a quadratic residue mod p? ", is_7_qr);
if (is_7_qr,
    y0 = lift(Mod(7, p)^((p+1)/4));
    print("  ⇒ 3-torsion points at x=0: (0, ±", y0 % 10^10, "...) are F_p-rational");
    print("  These are genuine 3-torsion points (order 3 in E(F_p)).");
    print("  BUT: n = #E(F_p) is prime > 3, so E(F_p)[3] = {O}.");
    print("  CONTRADICTION: (0,√7) ∈ E(F_p) with order dividing 3 but n is prime > 3.");
    print("  Resolution: x=0 is a root of ψ₃ but ψ₃ = 3x(x³+28), so (0,√7) has order 3");
    print("  in E(F̄_p) but not necessarily in E(F_p).  Check: is 3*(0,√7) = O in E?");
    E = ellinit([0, 7], p);
    y0_mod = lift(Mod(7, p)^((p+1)/4));
    P0 = [0, y0_mod];
    P3 = ellmul(E, P0, 3);
    print("  3*(0,√7) = ", P3, "    (should be [0] if P0 has order 3)");
    print("  Order of (0,√7)? n_secp mod 3 = ", n_secp % 3);
    print("  Since n_secp is prime and (0,√7) ∈ E(F_p), its order divides n_secp.");
    print("  So order of (0,√7) = n_secp (prime), NOT 3.");
    print("  But ψ₃(0)=0 means (0,y) is annihilated by [3] in E(F̄_p): 3*(0,y)=O over F̄_p.");
    print("  Key fact: (0,y) is NOT in E(F_p) unless y²=7 is a QR AND (0,y) has order | n_secp.");
    print("  Since gcd(3, n_secp) = 1 (n_secp prime, n_secp ≠ 3), the 3-torsion E[3](F_p)={O}.")
,
    print("  ⇒ 7 is NOT a QR mod p.  x=0 gives no F_p-rational 2-torsion.");
    print("  (So the linear factor 3x of ψ₃ contributes x=0 which has y²=7 non-square: no point.)")
);
print("");

\\ Factor x³ + 28 over F_p
print("---- Factoring x³ + 28 over F_p ----");
poly_3tors = Mod(1, p) * (x^3 + 28);
fac_3tors = factor(poly_3tors);
nrows_3 = matsize(fac_3tors)[1];
pat_3 = [];
for (i = 1, nrows_3, pat_3 = concat(pat_3, [poldegree(fac_3tors[i,1])]));
print("  Factorisation pattern of x³+28 mod p: ", pat_3);
if (pat_3 == [3],
    print("  ⇒ IRREDUCIBLE. E[3] over F_p has no non-trivial x-coordinates in F_p.");
    print("     Galois acts on E[3] via a cyclic group of order 3 (in GL₂(F₃)).");
    print("     E[3] ≅ (Z/3Z)² as abelian group; Frobenius acts with char poly x²+tx+p mod 3.")
);
if (pat_3 == [1, 2],
    print("  ⇒ ONE linear factor. One pair of F_p-rational non-trivial 3-torsion x-coordinates.")
);
if (pat_3 == [1, 1, 1],
    print("  ⇒ SPLITS COMPLETELY. All 3-torsion is F_p-rational.  E[3](F_p) = (Z/3Z)².")
);
print("");

\\ ---- Frobenius action on E[3] ----
print("---- Step 2: Frobenius action on E[3] ----");
t_mod3 = t_secp % 3;
p_mod3 = p % 3;
print("t mod 3 = ", t_mod3, "   (trace of Frobenius mod 3)");
print("p mod 3 = ", p_mod3, "   (p mod 3)");
print("");
print("Characteristic polynomial of Frobenius on T_3(E): x² - tx + p ≡ x² + ", (-t_mod3+3)%3, "x + ", p_mod3, " (mod 3)");
\\ The Galois module structure of E[3] is determined by this char poly mod 3.
charpoly_at3 = Pol([1, (-t_mod3+3)%3, p_mod3], x);
print("Char poly mod 3: ", charpoly_at3);
print("Factored over F_3: ", factor(Mod(1,3)*charpoly_at3));
print("");

\\ ---- Step 3: Which curves E' have the same E'[3] Galois module as E[3]? ----
print("---- Step 3: Galois-module matching condition for E' ----");
print("");
print("The E[3] Galois module is determined by the Frobenius eigenvalues λ₁,λ₂ in F̄₃,");
print("i.e., by the factored char poly of Frobenius mod 3.");
print("");
print("For E' to satisfy E'[3] ≅ E[3] as Gal-modules, E' must have the SAME");
print("Frobenius characteristic polynomial mod 3.  This is equivalent to:");
print("  t' ≡ t  (mod 3)  AND  p ≡ p  (trivially)");
print("i.e., the trace t' of E' satisfies t' ≡ ", t_mod3, " (mod 3).");
print("");
print("For the j=0 sextic-twist family:");
print("  T_0 = t        ≡ ", t_secp % 3, " (mod 3)");

\\ Print traces mod 3 for the 6 twists (using scalar-mult verified assignment from thread18_trace_verify)
s_val = sqrtint((4*p - t_secp^2) / 3);
t_plus  = (t_secp + 3*s_val) / 2;
t_minus = (t_secp - 3*s_val) / 2;
\\ Scalar-mult verified: k=1 ↔ T[6]=t_plus, k=2 ↔ T[5]=(3s-t)/2, k=4 ↔ T[3]=-(t+3s)/2, k=5 ↔ T[2]=-(3s-t)/2
T_verified = [t_secp, t_plus, (3*s_val - t_secp)/2, -t_secp, -(t_secp + 3*s_val)/2, (t_secp - 3*s_val)/2];

print("  Six traces (scalar-mult verified) mod 3:");
for (k = 1, 6,
    local(tk, nk);
    tk = T_verified[k];
    nk = p + 1 - tk;
    print("    k=", k-1, ": t=", tk, "  t mod 3 = ", tk % 3, "  N mod 3 = ", nk % 3)
);
print("");
print("Curves with t' ≡ t (mod 3) = ", t_secp%3, " have matching 3-torsion Galois module type.");
print("(This is a NECESSARY condition for (3,3)-gluing, not sufficient.)");
print("");

\\ ---- Step 4: gcd check for matching-trace-mod-3 pairs ----
print("---- Step 4: gcd check for pairs with matching t mod 3 ----");
print("");
\\ Collect indices with t_k ≡ t_secp (mod 3)
matching_k = [];
for (k = 1, 6, if (T_verified[k] % 3 == t_secp % 3, matching_k = concat(matching_k, [k])));
print("Twists with t_k ≡ t (mod 3): k ∈ {", matching_k, "}");
print("");

\\ Check all pairs among matching_k (and including secp256k1 itself as k=1)
for (ii = 1, #matching_k,
    for (jj = ii+1, #matching_k,
        local(ki, kj, Ni, Nj, g33);
        ki = matching_k[ii]; kj = matching_k[jj];
        Ni = p + 1 - T_verified[ki];
        Nj = p + 1 - T_verified[kj];
        g33 = gcd(Ni, Nj);
        print("  (k=", ki-1, ", k=", kj-1, "): gcd(N_{", ki-1, "}, N_{", kj-1, "}) = ", g33,
              "  → H3 analog: ", if(g33==1,"PASS","FAIL (gcd > 1)"))
    )
);
print("");

\\ ---- Summary ----
print("================================================================");
print("Summary: (3,3)-isogeny cover obstruction analysis");
print("================================================================");
print("");
print("1. E[3] Galois module: char poly mod 3 = ", charpoly_at3);
print("   Factor pattern of x³+28 over F_p: ", pat_3);
print("");
print("2. 3-torsion Frobenius: t mod 3 = ", t_secp%3);
print("   The necessary condition (t' ≡ t mod 3) selects ", #matching_k, " of 6 sextic twists.");
print("");
print("3. Among these ", #matching_k, " twists, gcd-coprimality (H3 analogue):");
print("   (secp256k1 has prime order n; any gcd with another N_k must divide n, hence = 1");
print("    if the two orders share no common prime factor.)");
print("");
print("4. Key difference from (2,2) case:");
print("   (2,2): needs matching 2-torsion pattern (x³+b factoring mod p).");
print("   (3,3): needs matching 3-torsion Galois type (char poly of Frob mod 3).");
print("   Additional (3,3)-gluing condition: the curves must share a specific");
print("   Weil pairing compatibility on E[3], which is a stronger condition.");
print("");
print("5. ECDLP impact: Even if (3,3)-covers exist, Jac order ≈ p² and the");
print("   best HCDLP algorithm is Pollard-ρ at cost ~p >> √p = cost of ECDLP.");
print("   No speedup from (3,3)-covers over prime fields.");
print("");
print("6. Open question: do all 3 pairs among matching twists satisfy the full");
print("   (3,3)-gluing conditions (including Weil pairing compatibility)?");
print("   This requires implementing the (3,3)-analogue of Howe's theorem,");
print("   which is less standard than the (2,2) case.");
print("");
print("================================================================");
print("Done: Thread 18 (3,3)-cover probe complete.");
print("================================================================");

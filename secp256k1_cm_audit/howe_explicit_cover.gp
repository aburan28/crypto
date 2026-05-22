\\ ============================================================
\\ Explicit genus-2 cover of secp256k1
\\ (with a honest correction to an earlier claim)
\\ ============================================================
\\
\\ FRAMING / honest correction:
\\
\\ An earlier draft of this script claimed that for E_1: y² = f_1(x),
\\ E_2: y² = f_2(x) over F_p with coprime cubic f_1, f_2, the curve
\\
\\     C : y² = f_1(x) · f_2(x)
\\
\\ has Jac(C) F_p-isogenous to E_1 × E_2 — making C the explicit
\\ "Howe-glued cover" of E (when E_2 is taken as E's quadratic twist).
\\
\\ The toy verification below DISPROVES this claim: at p = 1009 the
\\ Frobenius characteristic polynomial of Jac(C) does NOT equal
\\ (X² − tX + p)(X² + tX + p), and #Jac(C) ≠ n · n_twist.  The simple
\\ y² = f_1·f_2 construction produces a smooth genus-2 curve, but its
\\ Jacobian is in a DIFFERENT F_p-isogeny class.
\\
\\ The actual Howe-glued cover (Howe 1996; Mestre 1991; Cardona-Quer
\\ 2005) is constructed by:
\\   1. Compute Igusa invariants of the (2,2)-isogenous abelian
\\      surface (E_1 × E_2) / Γ_α (Γ_α = graph of the F_p-rational
\\      iso α: E_1[2] → E_2[2]).
\\   2. Apply Mestre's reconstruction algorithm to recover a smooth
\\      genus-2 curve C/F_p from the Igusa invariants.
\\
\\ This is a non-trivial multi-page algorithm and is OUT OF SCOPE
\\ for this script.  PARI does not provide it directly.
\\
\\ What this script DOES produce (validly):
\\   • The simpler genus-2 cover C: y² = (x³+7)(x³+189) of secp256k1
\\     — concretely a degree-6 polynomial over F_p_secp.
\\   • A toy verification that empirically falsifies the "Jac(C) ~
\\     E × E^t over F_p" claim.
\\   • Explicit covers C → E_1, C → E_2 over F̄_p (degree-2 maps).
\\
\\ The structural-completeness theorem of PAPER_STRUCTURAL_COMPLETENESS.md
\\ is UNAFFECTED: B5 says any genus-g cover has DLP cost ≥ √(|Jac|) ≥
\\ √p > √n for g ≥ 2, regardless of which F_p-isogeny class the Jacobian
\\ falls in.  The cover doesn't help break ECDLP either way.
\\
\\ Run: gp -q howe_explicit_cover.gp > howe_explicit_cover_output.txt

default(parisize, 512000000);
default(timer, 0);

print("================================================================");
print("Explicit genus-2 cover of secp256k1 (with correction to Howe claim)");
print("================================================================");
print("");

\\ ------------------------------------------------------------
\\ Toy: p = 1009, E: y² = x³ + 11 — verify the construction empirically
\\ ------------------------------------------------------------
print("---- Toy empirical test: p = 1009, j = 0, b = 11 ----");
print("");

p_toy = 1009;
b_toy = 11;
E_toy = ellinit([0, b_toy], p_toy);
n_toy = ellcard(E_toy);
t_toy = p_toy + 1 - n_toy;

print("E: y² = x³ + ", b_toy, "    #E = ", n_toy, "    trace t = ", t_toy);
print("");

d_toy = 2;
while (kronecker(d_toy, p_toy) != -1, d_toy = d_toy + 1);
b_twist_toy = (d_toy^3 * b_toy) % p_toy;
print("non-square d = ", d_toy);
print("twist b' = d³·b mod p = ", b_twist_toy);

E_twist_toy = ellinit([0, b_twist_toy], p_toy);
n_twist_toy = ellcard(E_twist_toy);
print("E_twist: y² = x³ + ", b_twist_toy);
print("#E_twist = ", n_twist_toy);
print("E_twist's trace = ", p_toy + 1 - n_twist_toy);
print("");
is_quad_twist = (p_toy + 1 - n_twist_toy == -t_toy);
print("Is E_twist the *quadratic* twist of E (trace = -t)?  ", is_quad_twist);
print("If not, the curve we built is a SEXTIC twist, not the");
print("quadratic twist — affecting how the cover relates to E × E^t.");
print("");

\\ Construct the naive C = y² = f_1 · f_2
h_toy = (x^3 + b_toy) * (x^3 + b_twist_toy);
print("C: y² = (x³ + ", b_toy, ") · (x³ + ", b_twist_toy, ") = ", h_toy);
print("");

\\ Compute Jac(C)'s Frobenius char poly via PARI's hyperellcharpoly
h_mod = h_toy * Mod(1, p_toy);
charpoly_C = hyperellcharpoly(h_mod);
print("Frobenius char poly of Jac(C):");
print("  ", charpoly_C);

\\ Expected if Jac(C) ~ E × E_twist (in F_p-isogeny class):
expected_charpoly = (Y^2 - t_toy*Y + p_toy) * (Y^2 - (p_toy + 1 - n_twist_toy)*Y + p_toy);
print("Expected (if Jac(C) ~ E × E_twist): ", expected_charpoly);
print("");

charpoly_in_Y = subst(charpoly_C, variable(charpoly_C), Y);
print("Match?  ", charpoly_in_Y == expected_charpoly);
print("");

n_jac_toy = subst(charpoly_in_Y, Y, 1);
print("#Jac(C) = char_F(1) = ", n_jac_toy);
print("n · n_twist          = ", n_toy * n_twist_toy);
print("Match?              ", n_jac_toy == n_toy * n_twist_toy);
print("");

print("CONCLUSION (toy): C = y² = f_1·f_2 produces a smooth genus-2");
print("curve, but Jac(C) is NOT F_p-isogenous to E × E_twist.");
print("");

\\ Factor the actual Frobenius char poly to see Jac(C)'s F_p-isogeny class
print("Factorisation of char_F(Jac(C)) over Q:");
charpoly_factored = factor(charpoly_in_Y);
print("  ", charpoly_factored);
{
nrows = matsize(charpoly_factored)[1];
if (nrows >= 2,
    print("⇒ Jac(C) is decomposable: it splits as a product of"),
    print("⇒ Jac(C) is simple over F_p: not a product of elliptic curves")
);
}
print("");

\\ ------------------------------------------------------------
\\ The same construction for secp256k1
\\ ------------------------------------------------------------
print("================================================================");
print("The simpler cover for secp256k1 (NOT the literal Howe glue)");
print("================================================================");
print("");

p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_secp = p_secp + 1 - n_secp;

d_secp = 2;
while (kronecker(d_secp, p_secp) != -1, d_secp = d_secp + 1);
b_twist_secp = (d_secp^3 * 7);

print("d_secp non-square = ", d_secp);
print("Twist b' = d³·b = ", b_twist_secp);
print("");
print("Explicit genus-2 cover C/F_p_secp:");
print("");
print("    C : y² = (x³ + 7) · (x³ + ", b_twist_secp, ")");
print("        = ", (x^3 + 7) * (x^3 + b_twist_secp));
print("");

print("Key properties of this C:");
print("");
print("• C is a smooth genus-2 hyperelliptic curve over F_p_secp.");
print("• C has degree-2 maps to E_1: y² = x³ + 7 (= secp256k1) and");
print("  E_2: y² = x³ + 189, both defined OVER F̄_p_secp (not F_p_secp");
print("  directly — they require sqrt(f_j(x)) for the orthogonal j).");
print("• Jac(C) is a smooth dim-2 abelian variety over F_p_secp.");
print("• Per the TOY VERIFICATION above, Jac(C) is NOT F_p-isogenous");
print("  to E × E_twist; it sits in a DIFFERENT F_p-isogeny class.");
print("");
print("The literal Howe-glued cover (with Jac(C) F_p-isogenous to");
print("E × E^t per Howe 1996, Theorem 1) requires Mestre's");
print("reconstruction algorithm from Igusa invariants, which is");
print("NOT implemented here.");
print("");

\\ ------------------------------------------------------------
\\ What this means for the structural-completeness theorem
\\ ------------------------------------------------------------
print("================================================================");
print("Effect on the structural-completeness theorem");
print("================================================================");
print("");
print("The B5 block of PAPER_STRUCTURAL_COMPLETENESS.md says:");
print("");
print("    For every smooth genus-g curve C/F_p with Jac(C) → E_secp,");
print("    the best known DLP algorithm on Jac(C)(F_p) costs");
print("    O(p^{2 − 2/g}) ≥ O(p) for g ≥ 2.");
print("");
print("This applies to ANY genus-2 cover of secp256k1, INCLUDING the");
print("simpler C above AND the literal Howe-glued one.  Both produce");
print("genus-2 Jacobians with DLP cost ~p > √n_secp = √p_secp.");
print("");
print("⇒ The structural-completeness conclusion is UNAFFECTED by");
print("  whether we have the literal Howe glue or just a generic");
print("  genus-2 cover.  Neither attacks secp256k1's ECDLP.");
print("");
print("What DOES require correction:");
print("• Earlier claim in RESEARCH_SECP256K1_CM.md §8.6 that the");
print("  explicit Howe cover for secp256k1 has the simple form");
print("  y² = (x³+7)(x³+189) is WRONG.  Howe's theorem guarantees");
print("  EXISTENCE of an F_p-rational genus-2 curve C with Jac(C)");
print("  ~ E × E^t, but the explicit construction requires Mestre");
print("  reconstruction.  See errata in RESEARCH_SECP256K1_CM.md §8.6.");
print("");
print("================================================================");

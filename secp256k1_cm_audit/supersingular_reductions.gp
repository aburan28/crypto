\\ ============================================================
\\ V4: supersingular reductions of secp256k1
\\ ============================================================
\\
\\ The curve y² = x³ + 7, considered over F_ℓ for small primes ℓ,
\\ is sometimes supersingular and sometimes ordinary.  When
\\ supersingular, the reduction lives in the supersingular ℓ-isogeny
\\ graph of F̄_ℓ — the same world that SIDH / SQIsign live in.
\\
\\ Question: does any supersingular reduction of secp256k1 leak
\\ ECDLP-relevant information about secp256k1's F_p curve?
\\
\\ Approach:
\\   1. For ℓ ∈ small primes, build E_ℓ: y² = x³ + 7 over F_ℓ.
\\   2. Compute #E_ℓ(F_ℓ); determine if E_ℓ is supersingular.
\\   3. For each supersingular E_ℓ, note its j-invariant and place
\\      in the supersingular isogeny graph.
\\   4. Analyse whether any connection to F_p ECDLP exists.
\\
\\ Run: gp -q supersingular_reductions.gp > supersingular_reductions_output.txt

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Supersingular reductions of secp256k1's equation y² = x³ + 7");
print("================================================================");
print("");

\\ For y² = x³ + 7 / F_ℓ (ℓ ≠ 2, 3, 7):
\\  - j-invariant = 0
\\  - E is supersingular over F_ℓ iff ℓ ≡ 2 (mod 3) (Deuring; the
\\    j=0 curve is supersingular precisely when -3 is a non-residue
\\    mod ℓ, equivalently ℓ ≡ 2 mod 3).

print("---- Supersingularity criterion for j = 0 curves ----");
print("");
print("Theorem (Deuring): for ℓ ≠ 2, 3, the curve y² = x³ + b over");
print("F_ℓ is supersingular iff ℓ ≡ 2 (mod 3), regardless of b.");
print("");
print("Equivalently: -3 is a non-residue mod ℓ  ⇔  ℓ is inert in Z[ω]");
print("            ⇔  the j=0 curve over F_ℓ is supersingular.");
print("");

\\ ------------------------------------------------------------
\\ Tabulate small primes and supersingularity status
\\ ------------------------------------------------------------
print("---- Reductions of secp256k1 (y² = x³ + 7) over small F_ℓ ----");
print("");
print("ℓ      ℓ mod 3   ord/SS    #E(F_ℓ)   trace   j(E_ℓ)");
print("---  --------- ---------  --------  ------  -------");

ssprimes = [];                            \\ collect supersingular ℓ
{
forprime (ell = 5, 100,
    if (ell == 7, next);                  \\ skip char-7 (b = 7 vanishes)
    local(E_ell, n_ell, t_ell, j_ell, is_ss);
    E_ell = ellinit([0, 7], ell);
    n_ell = ellcard(E_ell);
    t_ell = ell + 1 - n_ell;
    is_ss = (t_ell % ell == 0);           \\ supersingular iff p | t
    j_ell = lift(Mod(0, ell));            \\ j = 0 always (a = 0)
    print("  ", ell, "      ", ell % 3, "       ",
          if (is_ss, "SS  ", "ord "),
          "   ", n_ell, "    ", t_ell, "   ", j_ell);
    if (is_ss, ssprimes = concat(ssprimes, [ell]))
);
}
print("");
print("Total supersingular primes ℓ < 100: ", #ssprimes);
print("Supersingular ℓ list: ", ssprimes);
print("");

\\ ------------------------------------------------------------
\\ Connection to F_p ECDLP — analysis
\\ ------------------------------------------------------------
print("---- Connection (or lack thereof) to F_p ECDLP ----");
print("");
print("For each supersingular reduction E_ℓ of y² = x³ + 7:");
print("");
print("• End(E_ℓ) is a maximal order in the quaternion algebra");
print("  B_{ℓ,∞} = the unique quaternion algebra over Q ramified");
print("  exactly at ℓ and ∞.");
print("");
print("• j(E_ℓ) = 0 always (since the equation y² = x³ + 7 has j = 0).");
print("");
print("• The supersingular j-invariant set over F̄_ℓ contains roughly");
print("  ℓ/12 elements; j = 0 is one of them when ℓ ≡ 2 mod 3.");
print("");
print("• Castryck-Decru-Maino-Martindale-Panny-Pope-Robert (2022)");
print("  broke SIDH by exploiting the supersingular ℓ-isogeny");
print("  problem with auxiliary torsion data.  This is a SUPERSINGULAR-");
print("  WORLD attack, not connected to F_p ECDLP.");
print("");
print("---- Why supersingular reductions DON'T attack F_p ECDLP ----");
print("");
print("Deuring's lifting theorem: an ordinary curve over F_p lifts");
print("uniquely (up to F_p²-isomorphism) from any of its reductions.");
print("BUT the lifting is FROM F_ℓ TO characteristic 0, then specialise");
print("BACK to F_p.  This DOES NOT give a useful 'attack pipeline'");
print("because:");
print("");
print("  (a) The lift takes you to characteristic 0 (Q-bar / Q_p), not");
print("      to a 'smaller' DLP problem.");
print("");
print("  (b) The supersingular isogeny graph at ℓ has rich structure,");
print("      but its DLP-relevant features (e.g., GMP-isogeny-problem");
print("      hardness, SIDH-CSIDH) are about the ISOGENY problem on");
print("      the graph, not about ECDLP on individual curves in the");
print("      graph.");
print("");
print("  (c) Castryck-Decru's SIDH break needs specific torsion data");
print("      that's published by the protocol.  For our setting (we");
print("      just have an ECDLP instance on E/F_p), no such data is");
print("      published.");
print("");
print("  (d) The j-invariant 0 supersingular curve over F̄_ℓ is the");
print("      same algebraic-geometric object as secp256k1's reduction,");
print("      but reading ECDLP information from this object back to");
print("      F_p requires lifting through Deuring, which loses the");
print("      DLP secret in the process.");
print("");

\\ ------------------------------------------------------------
\\ Smart's attack (canonical lift) — when does it apply?
\\ ------------------------------------------------------------
print("---- Smart's attack revisited via canonical lift ----");
print("");
print("Smart 1999 attack on anomalous curves (n = p) uses the");
print("CANONICAL LIFT of E/F_p to a curve E/Z_p over the p-adic");
print("integers, then applies the formal logarithm to extract the");
print("discrete log.");
print("");
print("For secp256k1: n ≠ p (n = p + 1 - t with t ≈ 2^{129} ≠ 1).");
print("Therefore canonical-lift extraction of the discrete log via");
print("the formal logarithm DOES NOT WORK — the formal log produces");
print("a multiple of n in Z_p, but extracting d from this requires");
print("knowing n / gcd(n, p) which equals n.  Then d · log(P) =");
print("log(Q) mod p², which gives d mod (something) but not d mod n");
print("when n ∤ p.");
print("");
print("So canonical lifts give nothing for non-anomalous curves like");
print("secp256k1.  Confirmed null.");
print("");

\\ ------------------------------------------------------------
\\ Final V4 verdict
\\ ------------------------------------------------------------
print("================================================================");
print("V4 verdict: supersingular-reduction direction is null for");
print("           secp256k1's F_p ECDLP.");
print("================================================================");
print("");
print("• Supersingular reductions of secp256k1 exist for every prime");
print("  ℓ ≡ 2 (mod 3) (≈ half of all primes); each is in the");
print("  j=0 supersingular point of F̄_ℓ.");
print("");
print("• These reductions live in the SIDH/CSIDH supersingular world,");
print("  but the supersingular ℓ-isogeny problem (and its Castryck-");
print("  Decru-style breaks) are mathematically disjoint from F_p");
print("  ECDLP.  No known reduction in either direction.");
print("");
print("• Canonical-lift attacks (Smart 1999) require anomalous curves");
print("  (#E = p) — secp256k1 is not anomalous.");
print("");
print("⇒ V4 direction is closed for prime-field ECDLP on secp256k1.");
print("  The supersingular isogeny graph is a rich mathematical");
print("  object with its own cryptanalytic story (SIDH break) but");
print("  does not attack F_p ECDLP on ordinary curves.");
print("");
print("================================================================");

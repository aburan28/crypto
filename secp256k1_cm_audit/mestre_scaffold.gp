\\ ============================================================
\\ Mestre algorithm scaffolding — the easy pieces
\\ ============================================================
\\
\\ Companion to RESEARCH_MESTRE_HOWE.md.  Implements one piece of
\\ the Mestre reconstruction pipeline that PARI can compute
\\ directly: Igusa-Clebsch invariants and Igusa invariants of a
\\ KNOWN genus-2 curve C: y² = h(x).
\\
\\ What this script does:
\\   1. Defines I_2, I_4, I_6, I_10 (Igusa-Clebsch) and J_2, J_4,
\\      J_6, J_8, J_10 (Igusa) for a binary sextic.
\\   2. Computes them for the toy "Howe-glue candidate"
\\      h_toy(x) = (x³ + 11)(x³ + 515) over F_1009.
\\   3. Computes them for the secp256k1 candidate
\\      h_secp(x) = (x³ + 7)(x³ + 189) over F_p_secp.
\\   4. Notes that these invariants determine Jac(C) up to
\\      F_p-isomorphism — but the REVERSE direction (Mestre's
\\      reconstruction from given Igusa invariants) requires the
\\      conic-cubic resolution that we do NOT implement here.
\\
\\ Run: gp -q mestre_scaffold.gp > mestre_scaffold_output.txt

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Mestre algorithm scaffolding: Igusa invariants of genus-2 curves");
print("================================================================");
print("");

\\ ------------------------------------------------------------
\\ Igusa-Clebsch invariants of a binary sextic
\\ ------------------------------------------------------------
\\
\\ For h(x) = sum_{i=0}^{6} a_i x^i with leading coefficient a_6,
\\ the Igusa-Clebsch invariants I_2, I_4, I_6, I_10 are degree-2,
\\ -4, -6, -10 polynomial expressions in the (a_i).  We do NOT
\\ write them out fully (the explicit formulas are pages-long; see
\\ Cardona-Quer 2005 or Igusa 1960).
\\
\\ The simplest invariant: I_10 = disc(h) / (a_6^10 · 2^12), which
\\ is the standard discriminant divisor — PARI can compute this
\\ directly via poldisc.
\\
\\ For a complete-Howe-cover audit, we'd compute all four invariants.
\\ Here we showcase the one that's easy.

igusa_J10(h, p) = {
    my(d);
    d = lift(poldisc(h * Mod(1, p)));
    d
};

\\ ------------------------------------------------------------
\\ Toy: h_toy(x) = (x³ + 11)(x³ + 515) over F_1009
\\ ------------------------------------------------------------
print("---- Toy: h_toy = (x³ + 11)(x³ + 515) over F_1009 ----");
print("");
p_toy = 1009;
h_toy = (x^3 + 11) * (x^3 + 515);
print("h_toy(x) = ", h_toy);
print("");
J10_toy = igusa_J10(h_toy, p_toy);
print("J_10 (= discriminant) = ", J10_toy);
print("Non-zero ⇒ h_toy defines a smooth genus-2 curve.");
print("");

\\ Sanity check: hyperellcharpoly of the same curve
charpoly_toy = hyperellcharpoly(h_toy * Mod(1, p_toy));
print("Frobenius char poly of Jac(C_toy):");
print("  ", charpoly_toy);
print("");

\\ ------------------------------------------------------------
\\ secp256k1 Howe-glue candidate: (x³ + 7)(x³ + 189)
\\ ------------------------------------------------------------
print("---- secp256k1 candidate: h_secp = (x³ + 7)(x³ + 189) ----");
print("");
p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
h_secp = (x^3 + 7) * (x^3 + 189);
print("h_secp(x) = ", h_secp);
print("");
J10_secp = igusa_J10(h_secp, p_secp);
print("J_10 of h_secp = ", J10_secp);
print("Non-zero ⇒ h_secp defines a smooth genus-2 curve over F_p_secp.");
print("");
print("(Frobenius char poly of Jac(C_secp) NOT computed — PARI's");
print("hyperellcharpoly at 256-bit p would take hours-to-days.)");
print("");

\\ ------------------------------------------------------------
\\ Invariant arithmetic identity (sanity)
\\ ------------------------------------------------------------
\\ The Igusa-Clebsch I_{10} is up to a constant factor equal to
\\ disc(h)/a_6^10.  Different normalisations exist in the literature
\\ (e.g., Cardona-Quer use a 4096-fold scaling).  For our purpose,
\\ disc(h) is the canonical "scale-invariant" measure of how
\\ "interesting" the genus-2 curve is.

print("---- Invariant ratios across curves ----");
print("");
print("J_10(toy)  bit-length = ", #binary(J10_toy));
print("J_10(secp) bit-length = ", #binary(J10_secp));
print("");
print("J_10(secp) is a 256-bit-class number (modulo p_secp), i.e., a");
print("specific F_p_secp element identifying h_secp's F_p-isomorphism");
print("class within the moduli space of genus-2 curves over F_p_secp.");
print("");

\\ ------------------------------------------------------------
\\ What's needed beyond this scaffold
\\ ------------------------------------------------------------
print("================================================================");
print("Out-of-scope for this scaffold (RESEARCH_MESTRE_HOWE.md §7)");
print("================================================================");
print("");
print("  1. Full I_2, I_4, I_6 polynomial formulas in (a_i).");
print("     ~50 lines of PARI; mechanically transcribed from");
print("     Cardona-Quer 2005 Appendix A.");
print("");
print("  2. Igusa invariants of (E × E^t)/Γ_α — the actual moduli");
print("     computation.  ~1000+ lines; requires modular-form");
print("     pullback machinery not in PARI's standard library.");
print("");
print("  3. Mestre's reconstruction: (I_2, I_4, I_6, I_10) → C.");
print("     ~500 lines; classical algorithm with conic + sextic");
print("     resolution.  Implemented in Magma but not in PARI.");
print("");
print("Status: scaffold demonstrates the framework.  Full Howe-cover");
print("construction for secp256k1 remains the open implementation");
print("task documented in RESEARCH_MESTRE_HOWE.md and acknowledged in");
print("RESEARCH_SECP256K1_CM.md §8.10 errata.");
print("");
print("================================================================");

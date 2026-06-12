\\ igusa_f43_howe.gp
\\ Full Igusa-Clebsch tuple (J2,J4,J6,J10) for the Howe-glued curve
\\ C: y^2 = x^6 + 20x^3 + 5 over F_43.
\\
\\ 2026-06-11 computed only J2=41 and J10=36; this completes the tuple.
\\ Must be run from secp256k1_cm_audit/: gp -q igusa_f43_howe.gp

\\ Pre-set parisize to 256MB so igusa_clebsch_complete.gp's resize is a no-op
\\ (PARI 2.15.4: default(parisize,...) during read() aborts the read if it
\\  actually resizes the heap — documented in RESEARCH_AUTOLAB_LOG.md 2026-06-12)
default(parisize, 256000000);
read("igusa_clebsch_complete.gp");

p43 = 43;
h43 = x^6 + 20*x^3 + 5;

print("=== Part 2: Igusa (J2,J4,J6,J10) for C: y^2=x^6+20x^3+5 over F_43 ===");
print("(Howe-glued curve identified in 2026-06-11 session)");
print("");

igusa_res43 = igusa_quadruple(h43);

print("Igusa-Clebsch invariants over Z:");
print("  I2  = ", igusa_res43[1]);
print("  I4  = ", igusa_res43[2]);
print("  I6  = ", igusa_res43[3]);
print("  I10 = ", igusa_res43[4]);
print("");

I2v = igusa_res43[1] % p43;
I4v = igusa_res43[2] % p43;
I6v = igusa_res43[3] % p43;
I10v = igusa_res43[4] % p43;
if (I2v < 0, I2v = I2v + p43);
if (I4v < 0, I4v = I4v + p43);
if (I6v < 0, I6v = I6v + p43);
if (I10v < 0, I10v = I10v + p43);

print("Igusa-Clebsch mod 43:");
print("  J2  = ", I2v, "  (2026-06-11 reported J2=41:  match=", I2v == 41, ")");
print("  J4  = ", I4v, "  (NEW)");
print("  J6  = ", I6v, "  (NEW)");
print("  J10 = ", I10v, "  (2026-06-11 reported J10=36: match=", I10v == 36, ")");
print("");
print("Smooth (J10 != 0): ", I10v != 0);
print("");

\\ J8 from J2,J4,J6: J8=(J2*J6-J4^2)/4 (valid for char != 2)
J8v = lift(Mod(I2v*I6v - I4v^2, p43)) * lift(Mod(4,p43)^(-1));
J8v = lift(Mod(J8v, p43));
print("J8 = (J2*J6 - J4^2)/4 mod 43 = ", J8v);
print("");
print("Full Igusa tuple (J2:J4:J6:J8:J10) mod 43:");
print("  (", I2v, ":", I4v, ":", I6v, ":", J8v, ":", I10v, ")");
print("");
print("Cryptographic note: The Howe-glued curve y^2=x^6+20x^3+5");
print("over F_43 has Jac isogenous to E1 x E2 over F_{43^3}.");
print("The full Igusa tuple classifies it up to F_43-isomorphism.");

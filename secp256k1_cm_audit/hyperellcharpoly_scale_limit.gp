\\ hyperellcharpoly_scale_limit.gp
\\
\\ Documents a hard tooling limit discovered 2026-08-10: PARI's
\\ hyperellcharpoly (p-adic-Frobenius-based genus-2 point counting)
\\ overflows at the real secp256k1 prime (256 bits). This means ANY
\\ candidate Howe-cover verification against the real curve must go
\\ through algebraic Igusa-invariant matching (as chlrs_igusa_formula.gp
\\ already does), never naive point-counting -- there is no shortcut.
\\
\\ Run: gp -q hyperellcharpoly_scale_limit.gp

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
hh = Mod(1,p)*x^6 + Mod(196,p)*x^3 + Mod(1323,p);  \\ naive cover (x^3+7)(x^3+189)
print("Attempting hyperellcharpoly at the real secp256k1 prime (256 bits)...");
iferr(cp = hyperellcharpoly(hh), E, print("FAILED as expected: ", E));
print("");
print("Conclusion: hyperellcharpoly is infeasible above roughly toy-prime");
print("scale (worked fine at p=43, p=1009, p=1021 elsewhere in this dir).");
print("Point-counting-based verification of any Howe-cover formula candidate");
print("cannot be done directly on secp256k1; Igusa-invariant matching (see");
print("igusa_clebsch.gp, chlrs_igusa_formula.gp) is the only route at scale.");
quit;

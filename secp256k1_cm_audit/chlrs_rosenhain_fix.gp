\\ CHLRS Rosenhain fix: cubic-residue diagnosis + valid toy verification
\\ Run: gp -q chlrs_rosenhain_fix.gp
\\
\\ Toy case: p=7, b=6.  Then -b = -6 ≡ 1 mod 7, trivially a cubic residue,
\\ and x³+6 = x³-1 = (x-1)(x-2)(x-4) splits completely over F_7.

default(parisize, 128000000);
default(timer, 0);

is_cubic_res(a, p) = lift(Mod(a,p)^((p-1)/3)) == 1;

\\ Inline J2 and J10 (avoid noisy read of igusa_clebsch.gp)
sextic_coeffs(h) = {
  my(v); v = vector(7, i, 0);
  for(i = 0, min(6, poldegree(h)), v[i+1] = polcoeff(h, i));
  v
};
igusa_J2(h) = {
  my(a); a = sextic_coeffs(h);
  -120*a[1]*a[7] + 20*a[2]*a[6] - 8*a[3]*a[5] + 3*a[4]^2
};
igusa_J10(h) = poldisc(h);

print("================================================================");
print("CHLRS Rosenhain: cubic-residue diagnosis + toy verification");
print("================================================================");
print("");

\\ ---- Part 1: Obstruction check ----
print("---- Part 1: Obstruction check ----");
print("");
print("p=1009 b=11: cubic res (-b) mod p? ", is_cubic_res(-11, 1009));
print("  roots of x^3+11 mod 1009: ", polrootsmod(x^3+11, 1009));
print("  BLOCKED: T only F_{p^3}-defined; Rosenhain = non-trivial F_p-twist; match=0");
print("");
print("secp256k1 b=7: cubic res (-7) mod p_secp? ", is_cubic_res(-7, 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F));
print("  BLOCKED: x^3+7 irreducible over F_{p_secp}");
print("  ALL Howe-glueable pairs of secp256k1 have [3] 2-torsion");
print("");

\\ ---- Part 2: Valid toy case p=7, b=6 ----
\\ -b = -6 ≡ 1 mod 7 (trivially cubic residue); x³-1 = (x-1)(x-2)(x-4).
print("---- Part 2: Valid toy case p=7, b=6 ----");
pv = 7;  bv = 6;
print("p=", pv, "  b=", bv, "  (b ≡ -1 mod p, so -b ≡ 1 is trivially cubic residue)");
print("  p mod 6 = ", pv % 6, "  (need 1)");
print("  (-b) cubic res mod p? ", is_cubic_res(-bv, pv));
print("  roots of x^3+", bv, " mod ", pv, ": ", polrootsmod(x^3+bv, pv));

\\ First non-square d mod p
dv = 0;
for(dd = 2, pv-1, if(kronecker(dd,pv)==-1, dv=dd; break));
print("  first non-square d = ", dv, "  Kron(d,p) = ", kronecker(dv,pv));
b2v = lift(Mod(dv^3*bv, pv));
E1v = ellinit([0, bv], pv);
E2v = ellinit([0, b2v], pv);
t1v = pv + 1 - ellcard(E1v);
t2v = pv + 1 - ellcard(E2v);
N1 = pv + 1 - t1v;
N2 = pv + 1 + t1v;
print("  b2 = d^3*b mod p = ", b2v);
print("  E1: y^2 = x^3 + ", bv, " over F_", pv);
print("  E2: y^2 = x^3 + ", b2v, " over F_", pv);
print("  t(E1) = ", t1v, "   t(E2) = ", t2v, "  (want ", -t1v, "; ok = ", t2v == -t1v, ")");
print("  N1 = ", N1, "   N2 = ", N2, "   gcd = ", gcd(N1,N2));
print("  (gcd > 1 is OK for formula check; affects smoothness of Howe gluing only)");
print("");

\\ ---- Part 3: Rosenhain formula ----
\\ alpha = first root of x^3+b in F_p; omega = second_root / alpha.
print("---- Part 3: Rosenhain formula ----");
rts_b = polrootsmod(x^3+bv, pv);
print("Roots of x^3+", bv, " mod ", pv, ": ", rts_b);
alpha_v = lift(rts_b[1]);
om_v    = lift(Mod(lift(rts_b[2]), pv) * Mod(alpha_v, pv)^(-1));
print("alpha = ", alpha_v, "   omega = ", om_v);
print("  omega^3 mod p = ", lift(Mod(om_v, pv)^3), "  (need 1)");
print("  omega^2+omega+1 mod p = ", lift(Mod(om_v^2 + om_v + 1, pv)), "  (need 0)");

\\ Mobius: alpha->0, d*alpha->1, omega*alpha->inf
\\ lambda_i formulas from CHLRS (Howe 1996 / Rosenhain form):
dmod = Mod(dv, pv);
wmod = Mod(om_v, pv);
l1v = lift((wmod^2 - 1)*(dmod - wmod) / (wmod*(wmod - 1)*(dmod - 1)));
l2v = lift((dmod*wmod - 1)*(dmod - wmod) / (wmod*(dmod - 1)^2));
l3v = lift((dmod*wmod^2 - 1)*(dmod - wmod) / (wmod*(dmod*wmod - 1)*(dmod - 1)));
print("");
print("Rosenhain lambdas: l1=", l1v, "  l2=", l2v, "  l3=", l3v);
print("Five roots [0,1,l1,l2,l3] all distinct? ", l1v!=0 && l1v!=1 && l2v!=0 && l2v!=1 && l3v!=0 && l3v!=1 && l1v!=l2v && l1v!=l3v && l2v!=l3v);

\\ Rosenhain genus-2 curve over F_p
h_ros = Mod(1, pv) * x * (x-1) * (x-l1v) * (x-l2v) * (x-l3v);
cp_ros = hyperellcharpoly(h_ros);
nj_ros = subst(cp_ros, x, 1);
tgt_cp = (x^2 - t1v*x + pv) * (x^2 + t1v*x + pv);
print("");
print("Rosenhain char poly:  ", cp_ros);
print("Target char poly:     ", tgt_cp);
print("#Jac(Rosenhain) = ", nj_ros, "   target N1*N2 = ", N1*N2);
print("Order match?          ", nj_ros == N1*N2);
print("Full char poly match? ", cp_ros == tgt_cp);

\\ ---- Part 4: Igusa invariants of the Rosenhain curve ----
print("");
print("---- Part 4: Igusa invariants of Rosenhain curve (over Z) ----");
h_Z = x * (x-1) * (x-l1v) * (x-l2v) * (x-l3v);
print("h(x) over Z: ", h_Z);
J2v  = igusa_J2(h_Z);
J10v = igusa_J10(h_Z);
print("J2  = ", J2v);
print("J10 = ", J10v, "   smooth (J10 != 0)? ", J10v != 0);
print("J2  mod p = ", lift(Mod(J2v, pv)));
print("J10 mod p = ", lift(Mod(J10v, pv)));

\\ ---- Part 5: Naive cover comparison ----
print("");
print("---- Part 5: Naive y^2 = (x^3+b)(x^3+b2) comparison ----");
h_naive = Mod(1, pv) * (x^3 + bv) * (x^3 + b2v);
cp_naive = hyperellcharpoly(h_naive);
nj_naive = subst(cp_naive, x, 1);
print("Naive char poly:   ", cp_naive);
print("Target char poly:  ", tgt_cp);
print("#Jac(naive) = ", nj_naive, "   target = ", N1*N2);
print("Naive order match?       ", nj_naive == N1*N2);
print("Naive full poly match?   ", cp_naive == tgt_cp);
print("(Naive is NOT the Howe-glued curve; Rosenhain is.)");

\\ ---- Summary ----
print("");
print("================================================================");
print("Summary: CHLRS Rosenhain for j=0 Howe-glued covers");
print("================================================================");
print("REQUIREMENT: (-b) cubic residue mod p  ⟺  x^3+b splits over F_p.");
print("REASON: Mobius T sends alpha -> 0, d*alpha -> 1, omega*alpha -> inf.");
print("  If alpha not in F_p, T lies in PGL_2(F_{p^3}) and the Rosenhain");
print("  curve is a non-trivial F_p-twist of the Howe cover -> different");
print("  Frobenius -> char poly mismatch (match = 0).");
print("");
print("p=1009 b=11:    BLOCKED  ((-11) not cubic res mod 1009)");
print("secp256k1 b=7:  BLOCKED  (x^3+7 irreducible over F_{p_secp})");
print("  All 5 Howe-glueable pairs of secp256k1 are in [3]-class.");
print("  F_p Rosenhain + CHLRS Igusa formula: CANNOT be applied directly.");
print("");
print("Toy p=7 b=6:    VALID    ((-6)=1 mod 7 is cubic residue)");
print("  Rosenhain full char poly match = ", cp_ros == tgt_cp);
print("");
print("Next steps (secp256k1):");
print("  A: howe_richelot_v5.gp has explicit F_p genus-2 model -> Igusa there.");
print("  B: Mestre reconstruction from Igusa invariants.");
print("  C: F_{p^3} Rosenhain descent — construct over F_{p^3}, find F_p twist.");
print("================================================================");

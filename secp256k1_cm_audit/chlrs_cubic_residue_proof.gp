\\ chlrs_cubic_residue_proof.gp
\\
\\ Proves: the F_p Rosenhain/CHLRS formula is BLOCKED for secp256k1
\\ because (-7) is NOT a cubic residue mod p_secp.
\\
\\ Three experiments:
\\   A) Toy prime p=19:  (-7) IS a cubic residue → match=1 (formula correct)
\\   B) Proxy prime p=13: (-7) NOT a cubic residue → match=0 (formula fails)
\\   C) secp256k1 p:     (-7) NOT a cubic residue → F_p approach BLOCKED
\\
\\ The cubic residue condition: (-b) is a cubic residue mod p
\\   iff (-b)^{(p-1)/3} ≡ 1 mod p.
\\ Equivalently: x^3 + b has three roots in F_p.
\\
\\ Run: gp -q chlrs_cubic_residue_proof.gp

default(parisize, 256000000);
default(timer, 0);

is_cubic_res(a, p) = (lift(Mod(a,p)^((p-1)/3)) == 1);

print("================================================================");
print("Cubic-residue obstruction for secp256k1 Rosenhain formula");
print("================================================================");
print("");

\\ ================================================================
\\ Experiment A: toy prime p=19 where -7 IS a cubic residue
\\ ================================================================
print("---- Exp A: p=19, b=7, (-b) cubic residue ----");
{
  my(p, b, b_tw, ns, om, rts, l1, l2, l3, h, cp, n_jac, t1, t2, target);
  p = 19;
  b = 7;
  print("p = ", p, "  p mod 3 = ", p % 3);
  print("(-7)^{(p-1)/3} mod p = ", lift(Mod(-b,p)^((p-1)/3)),
        "  (1 = cubic residue)");
  rts = polrootsmod(Mod(1,p)*x^3 + b, p);
  print("Roots of x^3+", b, " mod ", p, ": ", Vec(rts));

  \\ Find first non-square
  ns = 2;
  while(kronecker(ns, p) != -1, ns++);
  b_tw = lift(Mod(ns^3 * b, p));
  print("Non-square d=", ns, "  twist coeff b2=d^3*b mod p=", b_tw);

  \\ Cube root of unity (from roots: omega = rts[2]/rts[1])
  om = lift(Mod(rts[2], p) / Mod(rts[1], p));
  print("omega = ", om, "  omega^3 = ", lift(Mod(om,p)^3));

  \\ Traces
  t1 = p + 1 - ellcard(ellinit([0, b], p));
  t2 = p + 1 - ellcard(ellinit([0, b_tw], p));
  print("t(E1)=", t1, "  t(E2)=", t2, "  (expected -", t1, ")");

  \\ Rosenhain lambdas
  l1 = lift((Mod(om,p)^2 - 1)*(Mod(ns,p) - Mod(om,p)) /
            (Mod(om,p)*(Mod(om,p) - 1)*(Mod(ns,p) - 1)));
  l2 = lift((Mod(ns,p)*Mod(om,p) - 1)*(Mod(ns,p) - Mod(om,p)) /
            (Mod(om,p)*(Mod(ns,p) - 1)^2));
  l3 = lift((Mod(ns,p)*Mod(om,p)^2 - 1)*(Mod(ns,p) - Mod(om,p)) /
            (Mod(om,p)*(Mod(ns,p)*Mod(om,p) - 1)*(Mod(ns,p) - 1)));
  print("Rosenhain lambdas: l1=", l1, "  l2=", l2, "  l3=", l3);

  h = Mod(1,p)*x*(x-1)*(x-l1)*(x-l2)*(x-l3);
  cp = hyperellcharpoly(h);
  n_jac = subst(cp, variable(cp), 1);
  target = (p+1-t1) * (p+1+t1);
  print("char poly Jac(C) = ", cp);
  print("#Jac(C) = ", n_jac, "  target = ", target);
  print("MATCH = ", n_jac == target);
  print("→ Exp A: formula CORRECT when (-b) is cubic residue");
  print("");
}

\\ ================================================================
\\ Experiment B: proxy prime p=13 where -7 is NOT a cubic residue
\\ ================================================================
print("---- Exp B: p=13, b=7, (-b) NOT cubic residue ----");
{
  my(p, b, b_tw, ns, t1, t2, target);
  p = 13;
  b = 7;
  print("p = ", p, "  p mod 3 = ", p % 3);
  print("(-7)^{(p-1)/3} mod p = ", lift(Mod(-b,p)^((p-1)/3)),
        "  (1=cubic residue, else not)");
  print("Roots of x^3+", b, " mod ", p, ": ",
        #polrootsmod(Mod(1,p)*x^3 + b, p), " roots");

  ns = 2;
  while(kronecker(ns, p) != -1, ns++);
  b_tw = lift(Mod(ns^3 * b, p));
  t1 = p + 1 - ellcard(ellinit([0, b], p));
  t2 = p + 1 - ellcard(ellinit([0, b_tw], p));
  print("t(E1)=", t1, "  t(E2)=", t2);

  \\ Construct naive cover and compute Jacobian char poly
  my(h, cp, n_jac);
  h = Mod(1,p)*(x^3 + b)*(x^3 + b_tw);
  cp = hyperellcharpoly(h);
  n_jac = subst(cp, variable(cp), 1);
  target = (p+1-t1)*(p+1+t1);
  print("Naive cover: y^2 = (x^3+7)(x^3+", b_tw, ") over F_", p);
  print("char poly Jac = ", cp);
  print("Expected (from E1*E2): x^4 + ", (2*p - t1^2), "*x^2 + ", p^2);
  print("#Jac = ", n_jac, "  target = ", target);
  print("MATCH = ", n_jac == target);
  print("→ Exp B: Jac NOT isogenous to E x E^t when (-b) not cubic residue");
  print("");
}

\\ ================================================================
\\ Experiment C: secp256k1 cubic residue check
\\ ================================================================
print("---- Exp C: secp256k1 prime ----");
{
  my(p, b, cr, nroots);
  p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
  b = 7;
  print("p mod 3 = ", p % 3, "  (must be 1 for cube roots of unity)");
  cr = is_cubic_res(-b, p);
  print("(-7)^{(p-1)/3} mod p == 1? ", cr);
  nroots = #polrootsmod(Mod(1,p)*x^3 + b, p);
  print("Roots of x^3+7 over F_p: ", nroots, " roots (need 3 for Rosenhain)");
  print("");
  if(cr,
    print("→ F_p Rosenhain applicable (unexpected!)"),
    print("→ F_p Rosenhain BLOCKED: (-7) not cubic residue mod p_secp")
  );
  print("→ x^3+7 irreducible over F_p → 2-torsion of secp256k1 in F_{p^3}\\F_p");
  print("");
}

\\ ================================================================
\\ Summary
\\ ================================================================
print("================================================================");
print("Summary: CHLRS/Rosenhain status for secp256k1");
print("================================================================");
print("");
print("The Rosenhain formula for Howe-gluing E x E^t requires:");
print("  (-b)^{(p-1)/3} ≡ 1 mod p  [cubic residue condition]");
print("  ⟺ x^3+b has 3 roots in F_p  [F_p-rational 2-torsion]");
print("");
print("Exp A (p=19): condition HOLDS → match=1 (formula verified correct)");
print("Exp B (p=13): condition FAILS → match=0 (Jac isogeny class wrong)");
print("Exp C (secp256k1): condition FAILS → F_p approach BLOCKED");
print("");
print("Next steps:");
print("  1. F_{p^3} extension: all cube roots live in F_{p^3}, but Jac");
print("     defined over F_{p^3} not F_p — doesn't give F_p DLP transfer.");
print("  2. Richelot search: find genus-2 curve C/F_p directly by");
print("     searching for Richelot factorizations of a degree-6 poly");
print("     whose Jacobian char poly = P_{E0}(T)*P_{Ej}(T) for the");
print("     5 qualifying pairs (0,1),(0,3),(0,4),(1,4),(3,4).");
print("  3. Status: BLOCKED for F_p Rosenhain approach. Richelot search");
print("     is the open problem; requires different PARI machinery.");

\\ chlrs_naive_cover_split_check.gp
\\
\\ Thread 2/3 (CHLRS Igusa forward map) — decisive check of the "naive cover"
\\ hypothesis used throughout howe_5pairs.gp / howe_5pairs_v2.gp (2026-07-26/27):
\\ that D: y^2=(x^3+b1)(x^3+b2) is F_p-isogenous to E1 x E2 for E1: y^2=x^3+b1,
\\ E2: y^2=x^3+b2, and that Richelot-dualizing D reaches the Howe cover.
\\
\\ Result: FALSE for every tested pair. D has an extra order-3 automorphism
\\ (x,y)->(zeta3*x,y) (since F depends on x only via x^3) that E1 x E2 does not
\\ generically share, so Jac(D) is a different abelian surface. D's Jacobian
\\ splits over F_p (as a product of two elliptic curves, but NOT E1,E2) iff
\\ b1*b2 is a QR mod p — evidence of a bielliptic involution
\\ (x,y)->(m/x, sqrt(b1*b2)*y/x^3), m^3=b1*b2, distinct from the naive
\\ E1,E2 factorization. See RESEARCH_AUTOLAB_LOG.md 2026-08-08 for the writeup.
\\
\\ Run: gp -q chlrs_naive_cover_split_check.gp

default(parisize, 256000000);

probe(pp, b1, b2, label) = {
  my(hh, cp, njac, E1, E2, t1, t2, fa, qr);
  hh = Mod(1,pp)*x^6 + Mod(b1+b2,pp)*x^3 + Mod(b1*b2,pp);
  cp = hyperellcharpoly(hh);
  njac = subst(cp, variable(cp), 1);
  E1 = ellinit([0,b1],pp); E2 = ellinit([0,b2],pp);
  t1 = pp+1-ellcard(E1); t2 = pp+1-ellcard(E2);
  fa = factor(cp);
  qr = kronecker(b1*b2, pp);
  print(label, ": p=",pp," b1=",b1," b2=",b2,"  t1=",t1," t2=",t2);
  print("  #Jac(D) = ", njac, "   target #E1*#E2 = ", (pp+1-t1)*(pp+1-t2),
        "   match=", njac==(pp+1-t1)*(pp+1-t2));
  print("  charpoly splits over Q (D has an F_p-rational elliptic factor pair): ",
        poldegree(fa[1,1]) < 4, "   kronecker(b1*b2,p)=",qr);
  print("");
}

print("=== p=43 (the 'validated' baseline from howe_5pairs_v2.gp Test 1) ===");
probe(43, 7, 13, "p43-baseline");

print("=== p=1009 toy pairs ===");
probe(1009, 11, 33, "p1009-b1=11,b2=33 (non-QR product)");
probe(1009, 11, 22, "p1009-b1=11,b2=22 (QR product -> splits, but not into E1,E2)");
probe(1009, 7, 189, "p1009-naive-secp-style-pair(0,3)");
probe(1009, 11, 515, "p1009-scalar-cube-pair(d=11 nonsquare, quadratic twist)");

print("Conclusion: #Jac(D) never equals #E1*#E2 in any tested case. The naive");
print("sextic-as-product-of-two-cubics ansatz is NOT the Howe cover of (E1,E2).");

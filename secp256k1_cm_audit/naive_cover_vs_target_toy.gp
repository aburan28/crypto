\\ naive_cover_vs_target_toy.gp
\\
\\ Directly checks (via hyperellcharpoly, not just Igusa invariants) whether
\\ the NAIVE cover y^2=(x^3+b1)(x^3+b2) has Jacobian ~ E1 x E2 for the
\\ non-degenerate secp256k1-type sextic-twist pairs (0,k), k=1..5, at a
\\ small toy prime. This check was not previously done directly against
\\ hyperellcharpoly for non-degenerate pairs (2026-07-26/27 sessions only
\\ checked the degenerate pair (0,3) via Igusa invariants).
\\
\\ PITFALL DISCOVERED (see RESEARCH_AUTOLAB_LOG.md 2026-08-10): p=1009,
\\ used throughout prior CHLRS toy-prime work, satisfies 36 | (p-1). This
\\ means the order-6 root of unity h (h^6=1) is CONTAINED in the 6th-power
\\ subgroup (F_p*)^6, so multiplying b1 by h^k for k=1..5 does NOT walk
\\ through 6 distinct sextic-twist classes at p=1009 -- it gives 5 curves
\\ all isomorphic to the original (same trace every time; verify below).
\\ Use p=1021 (36 does NOT divide 1020) for a valid toy-prime twist test.
\\
\\ Run: gp -q naive_cover_vs_target_toy.gp

check_pairs(pp, label) = {
  my(z3,h,b1);
  z3 = lift(polrootsmod(x^2+x+1,pp)[1]);
  h  = lift(polrootsmod(x^2-x+1,pp)[1]);
  b1 = 11;
  print("--- ", label, ": p=",pp,"  (p-1) mod 36 = ",(pp-1)%36," ---");
  for(k=1,5, b2 = lift(Mod(b1,pp)*Mod(h,pp)^k); E1 = ellinit([0,b1],pp); E2 = ellinit([0,b2],pp); t1 = pp+1-ellcard(E1); t2 = pp+1-ellcard(E2); target = (pp+1-t1)*(pp+1-t2); hh = Mod(1,pp)*x^6 + Mod(b1+b2,pp)*x^3 + Mod(b1*b2,pp); cp = hyperellcharpoly(hh); nj = subst(cp, variable(cp), 1); print("  pair (0,",k,"): b2=",b2,"  t1=",t1,"  t2=",t2,"  naive-Jac=",nj,"  target=",target,"  match=",nj==target));
};

check_pairs(1009, "DEGENERATE toy prime (36|p-1): twists collapse");
print("");
check_pairs(1021, "VALID toy prime: 6 distinct twist traces");
quit;

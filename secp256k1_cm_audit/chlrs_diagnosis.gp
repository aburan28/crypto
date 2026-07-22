\\ Diagnose CHLRS Rosenhain formula failure and find correct (2,2)-gluing
\\ for E1 x E2 -> Jac(C) on j=0 curves over F_p.

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("CHLRS diagnosis: finding correct Rosenhain model for Jac(C) ~ E1 x E2");
print("================================================================");
print("");

p = 1009; b = 11;
E1 = ellinit([0, b], p);
t1 = p + 1 - ellcard(E1);
print("E1: y^2 = x^3 + ", b, " over F_", p, "  trace = ", t1);

\\ Quadratic twist: E2 with trace -t1
\\ For j=0: quadratic twist of y^2=x^3+b is y^2=x^3+b*d^3 where d is non-square
d_ns = lift(Mod(2,p));
while(kronecker(d_ns,p) != -1, d_ns = d_ns+1);
b2 = lift(Mod(d_ns^3 * b, p));
E2 = ellinit([0, b2], p);
t2 = p + 1 - ellcard(E2);
print("E2: y^2 = x^3 + ", b2, " (d=", d_ns, ")  trace = ", t2);
print("Expected: t2 = ", -t1, " (quadratic twist)  OK? ", t2 == -t1);
print("");

\\ Target Weil polynomial for E1 x E2:
\\ (x^2 - t*x + p)(x^2 + t*x + p) = x^4 + (2p-t^2)*x^2 + p^2
t = t1;
target_mid = 2*p - t^2;
print("Target Weil poly: x^4 + ", target_mid, "*x^2 + ", p^2);
print("Target #Jac = ", 1 + target_mid + p^2);
print("");

\\ ---- Search: which y^2 = (x^3+b)(x^3+c) has the right char poly? ----
\\ Note: must search over F_p directly via Mod arithmetic.
print("---- Search over (x^3+b)(x^3+c), c=1..", p-1, " ----");
{
  local(c, h, cp, mid, ok);
  for(c = 1, p-1,
    \\ construct h over F_p
    h = Mod(1,p)*x^6 + Mod(c,p)*x^3 + Mod(b*c,p);
    if(poldisc(h) == 0, next);
    cp = hyperellcharpoly(h);
    \\ check if char poly is biquadratic (no x^3, x^1 terms)
    if(polcoeff(cp, 3) != 0 || polcoeff(cp, 1) != 0, next);
    mid = polcoeff(cp, 2);
    if(mid == Mod(target_mid, p),
      print("  MATCH: c=", c, "  char poly = ", cp);
    )
  );
}
print("Search done.");
print("");

\\ ---- Verify the Rosenhain formula from chlrs_igusa_formula.gp ----
\\ It gave char poly x^4 + 334*x^2 + 1018081 (wrong, expected 169)
\\ Let's check: 334 = 2*167, 169 = 13^2, 2*1009 - 13^2 = 2018 - 169 = 1849 = 43^2. Hmm.
\\ What curve has char poly x^4 + 334*x^2 + p^2?
\\ 334 = 2p - t'^2 => t'^2 = 2*1009 - 334 = 1684 = 4*421. sqrt(1684) is not integer.
\\ So the Rosenhain formula gives a SIMPLE Jacobian (not a product)!
print("---- What is 334 in x^4+334*x^2+p^2? ----");
print("334 = 2*", p, " - t^2 would need t^2 = ", 2*p-334, " = ", factor(2*p-334));
print("1684 is not a perfect square, so the Rosenhain Jac is SIMPLE (not a product)");
print("");

\\ ---- Try: which pairs (E, E_t) for j=0 over F_1009 give Howe-glueable products? ----
\\ The 6 j=0 curves over F_p have b in {b, b*g, b*g^2, b*g^3, b*g^4, b*g^5} (g = primitive root)
\\ where g is a primitive root mod p.
\\ From Thread 3: pair (0,3) [quadratic twist pair] is always glueable.
\\ Here (0,3) corresponds to b and b*g^3 = b*(-1)... wait, g^3 at order 6 is -g if g is prim root...
\\ Let me find all 6 sextic twists of y^2=x^3+1 for p=1009.
g_prim = znprimroot(p);
gen = lift(g_prim);
print("Primitive root mod ", p, " = ", gen);
bs = vector(6, k, lift(Mod(gen,p)^((p-1)/6 * (k-1))));
print("6 sextic-twist b-values: ", bs);
\\ For b=11, which index is it?
b_idx = -1;
for(k=1,6,
  \\ b=11 should be a sextic twist of b=1, i.e., b=11 = gen^((p-1)/6 * j) for some j?
  \\ Actually b=11 may not be in this list (it's a specific choice, not the canonical b=1)
  if(bs[k] == 11, b_idx = k)
);
print("b=11 is sextic twist index ", b_idx, " (1-based; -1 = not in canonical list)");
print("");

print("Script done. Key finding: -b=998 has (p-1)/3 power = ", lift(Mod(p-b,p)^((p-1)/3)),
      " (not 1, so -b is NOT a cube; 2-torsion of E1 not in F_p).");
print("This means the Rosenhain formula's cross-ratio construction does not map to E1 x E2 over F_p.");
print("Need: find which c gives char poly x^4 + ", target_mid, "*x^2 + p^2.");

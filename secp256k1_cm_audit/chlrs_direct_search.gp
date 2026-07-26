\\ chlrs_direct_search.gp
\\ Direct search for genus-2 curves with Jac (2,2)-isogenous to E1 x E2
\\ where E1, E2 are j=0 elliptic curves (quadratic twists of each other).
\\
\\ Run: gp -q chlrs_direct_search.gp

default(parisize, 256000000);
default(timer, 0);

p = 1009;
b1 = 11;

E1 = ellinit([0, b1], p);
t1 = ellap(E1, p);
p2 = p^2;
a2_target = 2*p - t1^2;

print("=============================================================");
print("Part 1: Direct search, p=1009, E1: y^2=x^3+11");
print("=============================================================");
print("E1 trace t = ", t1);
print("Target Weil poly: T^4 + ", a2_target, "*T^2 + ", p2);
print("Expected #Jac: ", 1 + a2_target + p2);
print("");

\\ ---- Helper: check char poly a2 coefficient ----
check_charpoly_a2(h) = {
  my(cp, c1, c2, c3);
  cp = hyperellcharpoly(h);
  c1 = polcoeff(cp, 3);
  c2 = polcoeff(cp, 2);
  c3 = polcoeff(cp, 1);
  if(c1 == 0 && c3 == 0, c2, -9999)
};

\\ ---- Part 1a: Sextic diagonal y^2 = x^6 + a*x^3 + b ----
print("---- Part 1a: Sextic diagonal y^2 = x^6 + a*x^3 + b ----");
print("Target a2 = ", a2_target);

search_sextic_diagonal() = {
  my(h, a2, found);
  found = List();
  forstep(a = 0, p-1, 1,
    forstep(b = 1, p-1, 1,
      h = Mod(1,p)*x^6 + Mod(a,p)*x^3 + Mod(b,p);
      a2 = check_charpoly_a2(h);
      if(a2 == a2_target, listput(found, [a, b]))
    )
  );
  found
};

found_sextic = search_sextic_diagonal();
print("Found ", #found_sextic, " sextic-diagonal curves with correct Weil poly.");
for(i=1, min(#found_sextic, 10),
  print("  a=", found_sextic[i][1], "  b=", found_sextic[i][2])
);
print("");

\\ ---- Part 1b: What does the Rosenhain Weil poly correspond to? ----
print("---- Part 1b: Diagnose Rosenhain output a2=334 ----");
\\ Rosenhain formula gave a2=334 instead of 169=a2_target.
\\ Find s1,s2 such that s1^2+s2^2 = 2p - 334 = 1684.
print("334 = 2p - s1^2 - s2^2 => s1^2 + s2^2 = ", 2*p - 334);
for(s = 0, 50,
  my(r);
  r = 1684 - s*s;
  if(r >= 0 && issquare(r), print("  s1=", s, " s2=", sqrtint(r)))
);
print("");

\\ ---- Part 1c: Find genus-2 curves with a2=334 (what Rosenhain gives) ----
print("---- Part 1c: Sextic diagonal curves with Rosenhain a2=334 ----");
for(a = 0, 10,
  for(b = 1, 50,
    my(h, a2v);
    h = Mod(1,p)*x^6 + Mod(a,p)*x^3 + Mod(b,p);
    a2v = check_charpoly_a2(h);
    if(a2v == 334, print("  a=", a, " b=", b))
  )
);
print("");

\\ ---- Part 2: Frobenius factorization of Rosenhain output ----
print("=============================================================");
print("Part 2: Factorize Weil poly T^4 + 334*T^2 + 1018081 over Q");
print("=============================================================");
\\ Does T^4 + 334*T^2 + p^2 factor as product of two quadratics?
\\ (T^2 + aT + p)(T^2 - aT + p) = T^4 + (2p-a^2)*T^2 + p^2
\\ a^2 = 2*1009 - 334 = 1684 (not a perfect square)
\\ So the Weil poly of the Rosenhain curve does NOT correspond to a product
\\ of two elliptic curves with the same trace ±a.
\\ But it might factor differently.

weil_rosenhain = x^4 + 334*x^2 + 1018081;
print("T^4 + 334*T^2 + p^2 factors over Q as:");
\\ Discriminant of x^2 + 334*x + 1018081 = 334^2 - 4*1018081 = 111556 - 4072324
disc = 334^2 - 4*p2;
print("  Discriminant (as quadratic in T^2) = 334^2 - 4*p^2 = ", disc);
\\ Roots T^2 = (-334 +/- sqrt(disc)) / 2
print("  T^2 = (-334 ± sqrt(", disc, ")) / 2");
if(issquare(-disc),
  print("  -disc is a perfect square: ", -disc, " = ", sqrtint(-disc), "^2"),
  print("  -disc = ", -disc, " (not a perfect square, Weil poly is irreducible over Q)")
);
print("");
\\ Over Q(sqrt(disc)): factor
sqd = sqrt(abs(disc)*1.0);
print("  Numerical: sqrt(-disc) = ", sqd, "i  =>  T^2 = (-334 ± ", sqd, "i) / 2");
\\ |T|^2 = p = 1009 => |T^2| = p, confirms these are Weil polynomial roots.
print("  |T^2| = sqrt(167^2 + (", sqd/2, ")^2) should = p=1009");
my(re, im);
re = -167;
im = sqd/2;
print("  |T^2| = ", sqrt(re^2 + im^2));
print("(Should be 1009 if these are valid Frobenius eigenvalues)");
print("");

\\ ---- Part 3: Summary ----
print("=============================================================");
print("Part 3: Summary");
print("=============================================================");

if(#found_sextic > 0, {
  my(a0, b0, cp0);
  a0 = found_sextic[1][1];
  b0 = found_sextic[1][2];
  cp0 = hyperellcharpoly(Mod(1,p)*x^6 + Mod(a0,p)*x^3 + Mod(b0,p));
  print("First correct curve: y^2 = x^6 + ", a0, "*x^3 + ", b0);
  print("Char poly: ", cp0);
  print("#Jac = ", subst(cp0, x, 1));
  print("Relationship to E1: y^2=x^3+", b1, " (trace=", t1, ")");
  print("  Note: if a0=0, b0=b1 then this is the 'naive cover' y^2=(x^3+b1)^2.");
  print("  If (x^3+b1)^2 factors: degenerate.");
  \\ Check if b0 = b1^2 mod p (naive factoring)
  if(b0 == lift(Mod(b1,p)^2),
    print("  b0 = b1^2 mod p: this is the squared form (degenerate)!"),
    print("  b0 != b1^2 mod p: genuine non-degenerate sextic.")
  );
  \\ Attempt to relate b0 to the twist structure
  print("  b1 = ", b1, "  b0/b1 mod p = ", lift(Mod(b0,p)*Mod(b1,p)^(-1)));
  print("  (b0/b1)^{(p-1)/2} = ", lift(Mod(b0*lift(Mod(b1,p)^(-1)),p)^((p-1)/2)),
        "  (1=square, -1=non-square)");
  \\ The b in E2 twist: b2 = d^3*b1 for some non-square d.
  \\ If the sextic is (x^3+b1)(x^3+b2) = x^6+(b1+b2)*x^3+b1*b2,
  \\ then a0 = b1+b2 mod p and b0 = b1*b2 mod p.
  a0_check = lift(Mod(b1 + lift(Mod(b0,p)*Mod(b1,p)^(-1)),p));
  print("  If sextic = (x^3+b1)(x^3+b2) then b2 = b0/b1 = ",
        lift(Mod(b0,p)*Mod(b1,p)^(-1)));
  my(b2_cand);
  b2_cand = lift(Mod(b0,p)*Mod(b1,p)^(-1));
  print("  Checking: a0 = b1+b2 mod p = ", lift(Mod(b1+b2_cand,p)));
  print("            a0 found = ", a0);
  \\ Compute trace of E2: y^2=x^3+b2
  my(E2_cand);
  E2_cand = ellinit([0, b2_cand], p);
  print("  Trace of E2 candidate (y^2=x^3+", b2_cand, "): ", ellap(E2_cand, p));
  print("  (For correct gluing we need trace = -t1 = ", -t1, ")")
}, {
  print("No sextic-diagonal curve found! The correct genus-2 curve is NOT")
  print("of the form y^2 = x^6 + a*x^3 + b.")
  print("The Howe-glued curve for j=0 pairs may require a more general form.")
});

print("");
print("=== DONE ===");

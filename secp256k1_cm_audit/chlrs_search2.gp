\\ chlrs_search2.gp — clean function-based search
\\ Finds genus-2 curves y^2 = x^6+a*x^3+b with Jac (2,2)-isogenous to E1xE2 (j=0 pair)

default(parisize, 256000000);
default(timer, 0);

p = 1009;
b1 = 11;
E1 = ellinit([0, b1], p);
t1 = ellap(E1, p);
a2_target = 2*p - t1^2;
print("p=",p," b1=",b1," t1=",t1," a2_target=",a2_target);

\\ Check if sextic y^2=x^6+a*x^3+b has the right char poly coefficient
ok_curve(a, b) = {
  my(h, cp, c1, c2);
  h = Mod(x^6 + a*x^3 + b, p);
  iferr(cp = hyperellcharpoly(h), e, return(0));
  c1 = polcoeff(cp, 3);
  c2 = polcoeff(cp, 2);
  (c1 == 0) && (c2 == a2_target)
};

\\ Return char poly for a curve
get_cp(a, b) = {
  my(h);
  h = Mod(x^6 + a*x^3 + b, p);
  hyperellcharpoly(h)
};

print("=== Searching sextic diagonal y^2 = x^6+a*x^3+b for a2=",a2_target," ===");
found_a = []; found_b = [];
for(a = 0, p-1,
  for(b = 1, p-1,
    if(ok_curve(a, b),
      found_a = concat(found_a, [a]);
      found_b = concat(found_b, [b])
    )
  )
);
print("Found ", #found_a, " curves.");
for(i = 1, min(#found_a, 15),
  my(ai, bi, cp);
  ai = found_a[i]; bi = found_b[i];
  cp = get_cp(ai, bi);
  print("  a=",ai,"  b=",bi,"  cp=",cp)
);
print("");

\\ Additional analysis: for each found curve, check if it has the split-sextic form
\\ y^2 = (x^3+b1)(x^3+b2) (i.e., a=b1+b2, b=b1*b2 mod p)
print("=== Checking split-sextic factorization (x^3+b1)(x^3+b2) ===");
for(i = 1, #found_a,
  my(ai, bi, b2_cand, t2);
  ai = found_a[i]; bi = found_b[i];
  \\ If split: ai = b1+b2 mod p, bi = b1*b2 mod p
  \\ => b2 = ai - b1 mod p (one candidate)
  b2_cand = lift(Mod(ai - b1, p));
  \\ Verify: b1*b2 == bi mod p
  if(lift(Mod(b1 * b2_cand, p)) == bi,
    my(E2c);
    E2c = ellinit([0, b2_cand], p);
    t2 = ellap(E2c, p);
    print("  Curve ",i,": SPLIT SEXTIC! b2=",b2_cand,
          "  trace(E2)=",t2," (want ",-t1,")")
  ,
    \\ Try other candidate: b2 = bi * b1^{-1}
    b2_cand = lift(Mod(bi, p) * Mod(b1, p)^(-1));
    if(lift(Mod(b1 + b2_cand, p)) == ai,
      my(E2c2);
      E2c2 = ellinit([0, b2_cand], p);
      t2 = ellap(E2c2, p);
      print("  Curve ",i,": SPLIT SEXTIC (b2=b/b1)! b2=",b2_cand,
            "  trace(E2)=",t2," (want ",-t1,")")
    ,
      print("  Curve ",i,": NOT split-sextic form (a=",ai,",b=",bi,")")
    )
  )
);
print("");

\\ Separately: check ALL non-square d and the corresponding (x^3+b1)(x^3+d^3*b1)
print("=== Split sextic y^2=(x^3+b1)(x^3+d^3*b1) for all non-square d ===");
found_d = [];
for(d = 2, p-1,
  if(kronecker(d, p) == -1,
    my(b2v, ai, bi, t2, E2v);
    b2v = lift(Mod(d, p)^3 * Mod(b1, p));
    ai = lift(Mod(b1 + b2v, p));
    bi = lift(Mod(b1 * b2v, p));
    if(ok_curve(ai, bi),
      E2v = ellinit([0, b2v], p);
      t2 = ellap(E2v, p);
      found_d = concat(found_d, [d]);
      print("  d=",d,"  b2=",b2v,"  trace(E2)=",t2,"  (want ",-t1,")")
    )
  )
);
print("Found ", #found_d, " valid non-square d values.");
print("");

\\ Report
print("=== Summary ===");
print("Target: Weil poly T^4+",a2_target,"*T^2+",p2," for p=",p);
print("#Found sextic-diagonal curves: ", #found_a);
if(#found_a > 0,
  print("First curve: y^2 = x^6+",found_a[1],"*x^3+",found_b[1])
);
print("Expected: the Howe-glued curve (E x E^t)/Gamma should be in this family.");
print("");
print("=== Weil poly diagnosis for Rosenhain formula output ===");
print("Rosenhain gave a2=334, our target is a2=",a2_target,"=",a2_target/169,"*169");
print("334 = 2*169 = 2*(2p-t^2). This is the a2 for E x E (same trace t,");
print("not the twist), OR for a different GL2(F2) kernel choice in the (2,2)-isogeny.");
print("");
print("DONE");

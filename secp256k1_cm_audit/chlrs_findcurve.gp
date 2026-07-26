\\ chlrs_findcurve.gp
\\ Find genus-2 curves (2,2)-isogenous to j=0 pairs via direct char poly search
\\ All multi-statement loops wrapped in functions

default(parisize, 256000000);
default(timer, 0);

p = 1009; b1 = 11;
E1 = ellinit([0, b1], p);
t1 = ellap(E1, p);
a2t = 2*p - t1^2;

print("p=",p,"  b1=",b1,"  t1=",t1,"  a2_target=",a2t);
print("Target Weil poly: T^4 + ",a2t,"*T^2 + ",p^2);

\\ Safely compute a2 coefficient of char poly of y^2 = f mod p
\\ f given as vector of coefficients [c0,c1,...,c5,c6] (constant to degree-6)
get_a2_sextic(a, b) = {
  my(h, cp);
  h = Mod(x^6 + a*x^3 + b, p);
  iferr(cp = hyperellcharpoly(h), e, return(-9999));
  if(polcoeff(cp, 3) != 0, return(-9999));
  polcoeff(cp, 2)
};

\\ Search ALL sextic-diagonal curves and print matches
search_sextic() = {
  my(cnt, a2v);
  cnt = 0;
  for(a = 0, p-1,
    for(b = 1, p-1,
      a2v = get_a2_sextic(a, b);
      if(a2v == a2t, cnt++; print("  a=",a,"  b=",b,"  a2=",a2v))
    )
  );
  print("Total: ",cnt)
};

\\ Only search split-sextic (x^3+b1)(x^3+b2) for b2 with trace(E2)=-t1
search_split() = {
  my(cnt, b2, E2c, t2, a_s, b_s, a2v);
  cnt = 0;
  for(b2 = 1, p-1,
    if(b2 != b1,
      E2c = ellinit([0, b2], p);
      t2 = ellap(E2c, p);
      if(t2 == -t1,
        a_s = lift(Mod(b1 + b2, p));
        b_s = lift(Mod(b1 * b2, p));
        a2v = get_a2_sextic(a_s, b_s);
        if(a2v == a2t,
          cnt++;
          print("  b2=",b2,"  a=",a_s,"  b=",b_s)
        )
      )
    )
  );
  print("Split-sextic count: ",cnt)
};

print("=== Split-sextic search (only twist pairs) ===");
search_split();
print("");

print("=== Full sextic-diagonal search (a in [0,p), b in [1,p)) ===");
print("NOTE: This takes ~30s for p=1009...");
search_sextic();
print("");
print("DONE");

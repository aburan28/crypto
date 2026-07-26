\\ chlrs_h3_analysis.gp
\\ Investigate Howe (H3) condition vs Rosenhain formula success.
\\ Run: gp -q chlrs_h3_analysis.gp

default(parisize, 256000000);
default(timer, 0);

mob(x_val, r1, r2, s1, p) =
  lift(Mod(x_val-r1,p) * Mod(x_val-r2,p)^(-1) * Mod(s1-r2,p) * Mod(s1-r1,p)^(-1));

test_row(p, b1, d_ns) = {
  my(b2, rts1, rts2, n1, n2, g, t1, r1, r2, r3, s1a, s2a, s3a, l1, l2, l3, h, cp, tgt);
  b2 = lift(Mod(b1*d_ns^3, p));
  rts1 = polrootsmod(Mod(1,p)*x^3 + b1, p);
  if (#rts1 < 3, return(-1));
  rts2 = polrootsmod(Mod(1,p)*x^3 + b2, p);
  if (#rts2 < 3, return(-1));
  n1 = ellcard(ellinit([0,b1], p));
  n2 = ellcard(ellinit([0,b2], p));
  t1 = p + 1 - n1;
  g = gcd(n1, n2);
  r1 = lift(rts1[1]); r2 = lift(rts1[2]); r3 = lift(rts1[3]);
  s1a = lift(rts2[1]); s2a = lift(rts2[2]); s3a = lift(rts2[3]);
  l1 = mob(r3, r1, r2, s1a, p);
  l2 = mob(s2a, r1, r2, s1a, p);
  l3 = mob(s3a, r1, r2, s1a, p);
  h = Mod(1,p)*x*(x-1)*(x-l1)*(x-l2)*(x-l3);
  cp = hyperellcharpoly(h);
  tgt = x^4 + (2*p - t1^2)*x^2 + p^2;
  [n1, n2, g, t1, cp==tgt]
};

print("p    b1  d   #E1   #E2   gcd  match H3?");

forprime(pp=7, 150,
  if(pp%3 != 1, next());
  for(bb=1, 8,
    my(dd, res);
    dd = 2;
    while(dd < pp && kronecker(dd,pp) != -1, dd++);
    if(dd >= pp, next());
    res = test_row(pp, bb, dd);
    if(res == -1, next());
    print(pp, "  ", bb, "  ", dd, "  ", res[1], "  ", res[2], "  ", res[3], "  ", res[5], "  ", if(res[3]==1,"H3:ok","H3:FAIL"))
  )
);

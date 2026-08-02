\\ naive_cover_falsification.gp
\\ ---------------------------------------------------------------------------
\\ Thread 2 companion — autolab 2026-08-02.
\\
\\ Several earlier entries treat the "naive cover"
\\        C_naive :  y^2 = (x^3 + b1 h^i)(x^3 + b1 h^j)
\\ as if it were the Howe-glued cover of the sextic-twist pair (E_i, E_j),
\\ and explain its failure on the pair (0,3) as a *degeneracy* (vanishing
\\ Richelot discriminant; 2026-07-27 entry, commit 9e9b0a8, Test 3).
\\
\\ This script tests the underlying premise directly, by point counting:
\\ does the Frobenius char poly of Jac(C_naive) equal that of E_i x E_j at all?
\\ (Char-poly equality, not merely #Jac = #E_i*#E_j, which is only necessary.)
\\
\\ Run:  gp -q naive_cover_falsification.gp
\\ ---------------------------------------------------------------------------

default(parisize, 500000000);
default(timer, 0);

\\ Frobenius char poly of y^2 = x^6 + aa*x^3 + bb, or 0 if singular.
cpjac(p, aa, bb) = {
  my(f);
  f = Mod(1,p)*x^6 + Mod(aa,p)*x^3 + Mod(bb,p);
  if(poldisc(f) == 0, return(0));
  subst(hyperellcharpoly(f), 'x, x);
}

\\ We test the CHAR-POLY identity, not just #Jac = #E_i*#E_j.  Order equality
\\ is only necessary, and on primes where the twist family collapses it is
\\ satisfied by accident (see the "distinct traces" column).
scan(p, b1) = {
  my(h, tr, ok, tot, A, B, aa, bb, cp, cp1, hits, ndist);
  h  = lift(znprimroot(p)^((p-1)/6));
  tr = vector(6, k, ellap(ellinit([0, lift(Mod(b1,p)*Mod(h,p)^(k-1))], p)));
  ndist = #Set(tr);
  ok = 0; tot = 0; hits = List();
  for(i = 0, 5,
    for(j = i+1, 5,
      A  = lift(Mod(b1,p)*Mod(h,p)^i);
      B  = lift(Mod(b1,p)*Mod(h,p)^j);
      aa = lift(Mod(A+B, p));
      bb = lift(Mod(A*B, p));
      cp = cpjac(p, aa, bb);
      if(cp == 0, next());
      cp1 = (x^2 - tr[i+1]*x + p) * (x^2 - tr[j+1]*x + p);
      tot++;
      if(cp == cp1, ok++; listput(hits, [i,j]))));
  print("  p=", p, "\tdistinct traces=", ndist,
        "\tnaive covers with correct char poly : ", ok, "/", tot,
        "\t", Vec(hits));
  [ok, tot, ndist];
}

print("===========================================================");
print(" Does the naive cover y^2=(x^3+b1 h^i)(x^3+b1 h^j) glue (E_i,E_j)?");
print("===========================================================");

TOY = [19,31,43,67,79,103,127,139,151,163,199,211,223,271,283,307,331,367,379,439,463,487,499];
{
  my(r, so, st, so6, st6, ndeg);
  so = 0; st = 0; so6 = 0; st6 = 0; ndeg = 0;
  for(k = 1, #TOY,
    r = scan(TOY[k], 7 % TOY[k]);
    so += r[1]; st += r[2];
    \\ "non-degenerate" = the 6 sextic twists really do have 6 distinct traces
    if(r[3] == 6, so6 += r[1]; st6 += r[2], ndeg++));
  print("");
  print("  AGGREGATE (all primes)          : ", so, " / ", st,
        "   = ", strprintf("%.1f", 100.0*so/st), " %");
  print("  AGGREGATE (6 distinct traces)   : ", so6, " / ", st6,
        "   = ", strprintf("%.1f", 100.0*so6/st6), " %");
  print("  (", ndeg, " primes dropped from the second line because the sextic");
  print("   twist family collapses to <6 trace values there, which makes");
  print("   accidental agreement common.)");
}

print("");
print("  Interpretation: the naive product cover is NOT the Howe cover of the");
print("  pair it is named after.  Inspecting the hit lists, essentially all of");
print("  them are the three quadratic-twist pairs (0,3),(1,4),(2,5) -- exactly");
print("  the a=0 cases y^2 = x^6 - A^2, where the theorem in");
print("  chlrs_forward_map.gp gives Jac ~ E_{-A^2} x E_{A^4}.  That coincides");
print("  with the named target pair only when -A^2 and A^4 happen to land in");
print("  the right sextic classes, which is a property of p, not a general");
print("  construction.");
print("");
print("  Consequence for the 2026-07-27 entry (commit 9e9b0a8, Test 3): the");
print("  failure on pair (0,3) was attributed to a Richelot degeneracy");
print("  (vanishing discriminant when d = h^3 = -1).  That diagnosis is");
print("  mis-scoped: the naive cover does not glue (E_0,E_3) on generic p in");
print("  the first place, so there is no degeneracy to repair.");

\\ ---------------------------------------------------------------------------
print("");
print("===========================================================");
print(" secp256k1: re-check the 2026-07-27 claim  189 = -7  mod 6th powers");
print("===========================================================");

P = 2^256 - 2^32 - 977;
{
  my(v, sq, cb);
  v  = lift(Mod(189,P) * Mod(-7,P)^(-1));      \\ = -27
  sq = (kronecker(v, P) == 1);
  cb = (lift(Mod(v,P)^((P-1)/3)) == 1);
  print("  189 / (-7) mod p = ", v);
  print("  is a square : ", sq);
  print("  is a cube   : ", cb);
  print("  is a 6th power : ", sq && cb);
  print("");
  if(sq && cb,
     print("  CONFIRMED: y^2=x^3+189 is F_p-isomorphic to y^2=x^3-7, so the"),
     print("  REFUTED: the 2026-07-27 identification does not hold."));
  if(sq && cb,
     print("  naive secp256k1 cover y^2=(x^3+7)(x^3+189) does correspond to"));
  if(sq && cb,
     print("  the index pair (0,3).  But per the toy-prime scan above, the"));
  if(sq && cb,
     print("  naive product cover does not glue (E_0,E_3) in the first place."));
}
print("");
print("done.");
quit

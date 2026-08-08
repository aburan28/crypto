\\ howe_forward_search.gp
\\
\\ Thread 3 (CHLRS Igusa forward map) — 2026-08-08 autolab.
\\
\\ 2026-07-27 (howe_5pairs_v2.gp Test 2) found that the "natural" 2-torsion
\\ identification beta = d*alpha (d = cube root of the twist ratio) gives the
\\ WRONG Jacobian order for p=1009, pair (E1: y^2=x^3+11, E2: y^2=x^3+515).
\\ The conclusion was that the forward map "requires the CHLRS Igusa inversion"
\\ i.e. is blocked without the literature formula.
\\
\\ This script tests a cheaper hypothesis first: the Richelot function is
\\ CORRECT (verified on the p=43 reference case) and the only bug is a WRONG
\\ CHOICE of Galois-equivariant identification alpha <-> beta. Since Frobenius
\\ acts on each 2-torsion Galois orbit {alpha, z3*alpha, z3^2*alpha} as a
\\ 3-cycle, the Galois-equivariant bijections to {beta, z3*beta, z3^2*beta}
\\ are exactly beta = d*z3^k*alpha for k=0,1,2 (3 choices, not 1). We only
\\ ever tried k=0. Brute-force try k=0,1,2 (and beta = -d*z3^k*alpha, since
\\ E2's OTHER independent x^3=rr' cube root set is {-beta} only in char != it's
\\ literally the same set — but the y-sign/Weil-pairing orientation gives an
\\ additional order-2 ambiguity, so also try swapping which of the two order-3
\\ covers of Frobenius we use) against the true target Jacobian order.
\\
\\ Run: gp -q howe_forward_search.gp

default(parisize, 256000000);

f3add(u, v, pp) = [(u[1]+v[1])%pp, (u[2]+v[2])%pp, (u[3]+v[3])%pp];
f3neg(u, pp)    = [(-u[1])%pp, (-u[2])%pp, (-u[3])%pp];
f3scl(c, u, pp) = [(c*u[1])%pp, (c*u[2])%pp, (c*u[3])%pp];

f3mul(u, v, rr, pp) = {
  my(c0, c1, c2, c3, c4);
  c0 = (u[1]*v[1]) % pp;
  c1 = (u[1]*v[2] + u[2]*v[1]) % pp;
  c2 = (u[1]*v[3] + u[2]*v[2] + u[3]*v[1]) % pp;
  c3 = (u[2]*v[3] + u[3]*v[2]) % pp;
  c4 = (u[3]*v[3]) % pp;
  [(c0 + rr*c3) % pp,  (c1 + rr*c4) % pp,  c2 % pp]
};

f3inv(u, rr, pp) = {
  my(a, b, c, nrm, ni);
  a = u[1]; b = u[2]; c = u[3];
  nrm = lift(Mod(a^3 + rr*b^3 + rr^2*c^3 - 3*rr*a*b*c, pp));
  if(nrm == 0, error("f3inv: norm=0, element not invertible"));
  ni = lift(Mod(nrm, pp)^(-1));
  [((a^2 - rr*b*c)*ni) % pp,
   ((rr*c^2 - a*b)*ni) % pp,
   ((b^2 - a*c)*ni) % pp]
};

richelot(sv, qv, z3, rr, pp) = {
  my(z3sq, G1c, G1x, G2c, G2x, G3c, G3x, D0, D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,lc,lcinv,aa,bb);

  z3sq = lift(Mod(z3, pp)^2);

  G1c = qv;          G1x = f3neg(sv, pp);
  G2c = f3scl(z3sq, qv, pp);  G2x = f3neg(f3scl(z3, sv, pp), pp);
  G3c = f3scl(z3, qv, pp);   G3x = f3neg(f3scl(z3sq, sv, pp), pp);

  D0 = f3add(f3add(
    f3add(f3mul(G2x,G3c,rr,pp), f3neg(f3mul(G3x,G2c,rr,pp),pp), pp),
    f3neg(f3mul(G1x, f3add(G3c, f3neg(G2c,pp), pp), rr, pp), pp), pp),
    f3mul(G1c, f3add(G3x, f3neg(G2x,pp), pp), rr, pp), pp);

  if(D0==[0,0,0], return([-1,-1]));
  D0inv = f3inv(D0, rr, pp);

  H1n2 = f3add(G3x, f3neg(G2x,pp), pp);
  H1n1 = f3scl(2, f3add(G3c, f3neg(G2c,pp), pp), pp);
  H1n0 = f3add(f3mul(G2x,G3c,rr,pp), f3neg(f3mul(G2c,G3x,rr,pp),pp), pp);
  H2n2 = f3add(G1x, f3neg(G3x,pp), pp);
  H2n1 = f3scl(2, f3add(G1c, f3neg(G3c,pp), pp), pp);
  H2n0 = f3add(f3mul(G3x,G1c,rr,pp), f3neg(f3mul(G3c,G1x,rr,pp),pp), pp);
  H3n2 = f3add(G2x, f3neg(G1x,pp), pp);
  H3n1 = f3scl(2, f3add(G2c, f3neg(G1c,pp), pp), pp);
  H3n0 = f3add(f3mul(G1x,G2c,rr,pp), f3neg(f3mul(G1c,G2x,rr,pp),pp), pp);

  H1x2=f3mul(H1n2,D0inv,rr,pp); H1x1=f3mul(H1n1,D0inv,rr,pp); H1x0=f3mul(H1n0,D0inv,rr,pp);
  H2x2=f3mul(H2n2,D0inv,rr,pp); H2x1=f3mul(H2n1,D0inv,rr,pp); H2x0=f3mul(H2n0,D0inv,rr,pp);
  H3x2=f3mul(H3n2,D0inv,rr,pp); H3x1=f3mul(H3n1,D0inv,rr,pp); H3x0=f3mul(H3n0,D0inv,rr,pp);

  P0=f3mul(H1x0,H2x0,rr,pp);
  P1=f3add(f3mul(H1x0,H2x1,rr,pp),f3mul(H1x1,H2x0,rr,pp),pp);
  P2=f3add(f3add(f3mul(H1x0,H2x2,rr,pp),f3mul(H1x1,H2x1,rr,pp),pp),f3mul(H1x2,H2x0,rr,pp),pp);
  P3=f3add(f3mul(H1x1,H2x2,rr,pp),f3mul(H1x2,H2x1,rr,pp),pp);
  P4=f3mul(H1x2,H2x2,rr,pp);

  Q0r=f3mul(P0,H3x0,rr,pp);
  Q1r=f3add(f3mul(P0,H3x1,rr,pp),f3mul(P1,H3x0,rr,pp),pp);
  Q2r=f3add(f3add(f3mul(P0,H3x2,rr,pp),f3mul(P1,H3x1,rr,pp),pp),f3mul(P2,H3x0,rr,pp),pp);
  Q3r=f3add(f3add(f3mul(P1,H3x2,rr,pp),f3mul(P2,H3x1,rr,pp),pp),f3mul(P3,H3x0,rr,pp),pp);
  Q4r=f3add(f3add(f3mul(P2,H3x2,rr,pp),f3mul(P3,H3x1,rr,pp),pp),f3mul(P4,H3x0,rr,pp),pp);
  Q5r=f3add(f3mul(P3,H3x2,rr,pp),f3mul(P4,H3x1,rr,pp),pp);
  Q6r=f3mul(P4,H3x2,rr,pp);

  if(Q0r[2]!=0||Q0r[3]!=0||Q3r[2]!=0||Q3r[3]!=0||Q6r[2]!=0||Q6r[3]!=0,return([-1,-1]));
  lc = Q6r[1]; if(lc==0, return([-1,-1]));
  lcinv = lift(Mod(lc,pp)^(-1));
  aa = (Q3r[1]*lcinv)%pp; bb=(Q0r[1]*lcinv)%pp;
  [aa, bb]
};

check_jac(aa, bb, target, pp) = {
  my(hh, cp, nj);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  [cp, nj, nj==target]
};

\\ ================================================================
\\ Try all Galois-equivariant identifications for a given (b1, b2, pp).
\\ rr = -b1 mod pp (so alpha^3 = rr are the x-coords of E1's 2-torsion).
\\ beta must satisfy beta^3 = -b2 mod pp; write beta = d*alpha for a FIXED
\\ cube root d of (b2/b1) mod pp (any one choice — the script already found
\\ one earlier), then the Galois-equivariant identifications are exactly
\\ beta_k = d*z3^k*alpha for k=0,1,2.
\\ ================================================================
try_all_gluings(b1, b2, t1_target, t2_target, pp, label) = {
  my(z3, d3, ratio, found, ok);
  print("--- ", label, "  (p=",pp,", b1=",b1,", b2=",b2,") ---");
  z3 = lift(polrootsmod(x^2+x+1,pp)[1]);
  ratio = lift(Mod(b2,pp)/Mod(b1,pp));
  \\ find a cube root of ratio in F_p (may not exist -> then d lives in F_{p^3},
  \\ skip, handled by a separate degree-9 search elsewhere)
  d3 = -1;
  for(cand=1, pp-1, if(lift(Mod(cand,pp)^3)==ratio, d3=cand; break));
  if(d3==-1,
    print("  ratio ", ratio, " has no cube root in F_p; skip (needs F_{p^3} alpha too)");
    return(0)
  );
  rr = (-b1)%pp;
  target = ((pp+1-t1_target)*(pp+1-t2_target)) % (2^300);
  target = (pp+1-t1_target)*(pp+1-t2_target);
  found = 0;
  for(k=0, 2,
    dk = lift(Mod(d3,pp)*Mod(z3,pp)^k);
    sv = [0, (1+dk)%pp, 0]; qv = [0, 0, dk%pp];
    res = richelot(sv, qv, z3, rr, pp);
    if(res[1]==-1,
      print("  k=",k," (beta=d*z3^",k,"*alpha): cover not defined over F_p"),
      chk = check_jac(res[1], res[2], target, pp);
      print("  k=",k," (beta=d*z3^",k,"*alpha): a=",res[1]," b=",res[2],
            "  #Jac=",chk[2],"  target=",target,"  match=",chk[3]);
      if(chk[3], found=1; print("    *** MATCH at k=",k," ***"))
    )
  );
  if(!found, print("  No k in {0,1,2} matches. Identification may need a sign flip on beta or an F_{p^3}-alpha too."));
  found
};

\\ ================================================================
print("================================================================");
print("Re-run Test 2 (p=1009) with all 3 Galois-equivariant gluings");
print("================================================================");
{
  pp=1009; b1=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  d=2; while(kronecker(d,pp)!=-1, d++);
  b2=lift(Mod(d^3*b1,pp));
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  E1: y^2=x^3+",b1,"  trace=",t1,"   E2: y^2=x^3+",b2,"  trace=",t2);
  try_all_gluings(b1, b2, t1, t2, pp, "Test2-p1009");
}
print("");

print("================================================================");
print("Sanity check: p=43 reference case (expect MATCH at some k)");
print("================================================================");
{
  pp=43; b1=7; b2=13;
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  E1: y^2=x^3+",b1,"  trace=",t1,"   E2: y^2=x^3+",b2,"  trace=",t2);
  try_all_gluings(b1, b2, t1, t2, pp, "p43-reference");
}
print("");

print("================================================================");
print("Scan: 10 fresh toy primes p=1 mod 3, p=3 mod 4 (irreducible x^3+b1),");
print("random non-square-cube d, check how often SOME k in {0,1,2} matches.");
print("================================================================");
{
  nprime=0; nmatch=0;
  forprime(pp=200, 2000,
    if(pp%3==1 && pp%4==3,
      b1=1;
      while(!issquarefree(x^3+b1)||polisirreducible(Mod(1,pp)*x^3+b1)==0, b1++;
        if(b1>pp,break));
      if(b1<=pp && polisirreducible(Mod(1,pp)*x^3+b1),
        d=2; while(kronecker(d,pp)!=-1 && d<pp, d++);
        if(d<pp,
          b2=lift(Mod(d^3*b1,pp));
          E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
          t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
          nprime++;
          r = try_all_gluings(b1, b2, t1, t2, pp, Str("scan p=",pp));
          if(r, nmatch++);
        )
      )
    )
  );
  print("");
  print("SCAN SUMMARY: ",nmatch,"/",nprime," primes had a matching k in {0,1,2}.");
}

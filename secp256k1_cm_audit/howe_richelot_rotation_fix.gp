\\ howe_richelot_rotation_fix.gp
\\
\\ Thread 3 (CHLRS Igusa forward map) — autolab 2026-08-09.
\\
\\ 2026-07-27's Test 2 (howe_5pairs_v2.gp) used the "naive" root pairing
\\ beta = d*alpha (a single fixed matching of E1's 3 non-trivial 2-torsion
\\ x-coords {alpha, z3*alpha, z3^2*alpha} to E2's {beta, z3*beta, z3^2*beta})
\\ and got the WRONG Jacobian order (1106283 vs target 1018251).
\\
\\ Hypothesis: a Galois-equivariant bijection E1[2]\{O} -> E2[2]\{O} is only
\\ constrained up to a *rotation* by the Z/3Z action (since Gal(F_{p^3}/F_p)
\\ permutes {alpha, z3*alpha, z3^2*alpha} cyclically and any Howe gluing map
\\ alpha_map must intertwine this action). There are exactly 3 such
\\ equivariant bijections: beta_i = z3^i * d * alpha, i=0,1,2. Test 2 only
\\ tried i=0. This script tries all three and checks which (if any) recovers
\\ the target Jacobian order (p+1-t1)*(p+1-t2).
\\
\\ Run: gp -q howe_richelot_rotation_fix.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---- F_{pp^3} arithmetic (copied from howe_5pairs_v2.gp) ----
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

check_jac(aa, bb, t_expected, pp) = {
  my(hh, cp, nj, target);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  target = (pp+1-t_expected) * (pp+1+t_expected);
  [cp, nj, target, nj==target]
};

print("================================================================");
print("p=1009, E1: y^2=x^3+11, E2: y^2=x^3+515 (d=11)  -- 3 rotations");
print("================================================================");
{
  pp=1009; b1=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  z3sq = lift(Mod(z3,pp)^2);
  d=2; while(kronecker(d,pp)!=-1, d++);
  b2=lift(Mod(d^3*b1,pp));
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  target = (pp+1-t1)*(pp+1-t2);
  print("  z3=",z3,"  d=",d,"  b2=",b2,"  t1=",t1,"  t2=",t2,"  target=#E1*#E2=",target);
  print("");
  rr=(-b1)%pp;
  rot = [1, z3, z3sq];
  for(i=1,3,
    di = lift(Mod(d*rot[i], pp));
    sv=[0,(1+di)%pp,0]; qv=[0,0,di%pp];
    res=richelot(sv,qv,z3,rr,pp);
    print("  rotation i=",i-1,":  beta=",di,"*alpha   Richelot -> a=",res[1],"  b=",res[2]);
    if(res[1]!=-1,
      chk=check_jac(res[1],res[2],t1,pp);
      print("    #Jac=",chk[2],"  target=",chk[3],"  match=",chk[4]);
      if(chk[4], print("    *** MATCH on rotation i=",i-1," ***"));
    ,
      print("    cover not defined over F_p for this rotation");
    );
  );
}
print("");
print("================================================================");
print("p=43 sanity re-check: does rotation i=0 still give the known-good");
print("Test-1 answer (a=41,b=5) when re-derived via the same d-based");
print("parameterization (b1=7,b2=13, so d^3=13/7 mod 43)?");
print("================================================================");
{
  pp=43; b1=7; b2=13;
  z3=6; z3sq=lift(Mod(z3,pp)^2);
  rr=(-b1)%pp;
  \\ find d with d^3 = b2/b1 mod pp
  dtarget = lift(Mod(b2,pp)/Mod(b1,pp));
  d=0; for(cand=1,pp-1, if(lift(Mod(cand,pp)^3)==dtarget, d=cand; break));
  print("  d (cube root of b2/b1=",dtarget,") = ",d);
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  target=(pp+1-t1)*(pp+1-t2);
  print("  t1=",t1,"  t2=",t2,"  target=",target);
  rot=[1,z3,z3sq];
  for(i=1,3,
    di=lift(Mod(d*rot[i],pp));
    sv=[0,(1+di)%pp,0]; qv=[0,0,di%pp];
    res=richelot(sv,qv,z3,rr,pp);
    print("  rotation i=",i-1,": a=",res[1],"  b=",res[2]);
    if(res[1]!=-1,
      chk=check_jac(res[1],res[2],t1,pp);
      print("    #Jac=",chk[2],"  target=",chk[3],"  match=",chk[4]);
    );
  );
}
print("");
print("Done.");

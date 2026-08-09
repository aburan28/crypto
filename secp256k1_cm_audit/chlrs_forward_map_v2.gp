\\ chlrs_forward_map_v2.gp
\\
\\ Thread 3 (CHLRS Igusa forward map), continued 2026-08-08.
\\ Literature access (arxiv.org, eprint.iacr.org, math.mit.edu) is EGRESS_BLOCKED
\\ in this environment -- cannot fetch Howe-Leprevost-Poonen or CHLRS directly.
\\ Instead: re-examine why the naive formula r1=-b1, r2=-b2 (validated at p=43,
\\ Test 1 of howe_5pairs_v2.gp) FAILED at p=1009 (Test 2, log 2026-07-27).
\\
\\ Hypothesis: Test 2 only tried ONE Galois embedding (sv=alpha+beta, qv=alpha*beta),
\\ but howe_richelot_v5.gp's do_class() tries THREE (beta, z3*beta, z3^2*beta) and
\\ only reports success if ANY of them lands in F_p. Test 2 may simply have picked
\\ the wrong branch. This script re-runs p=43 (sanity) and p=1009 with all 3 branches,
\\ AND brute-forces over the 3x3 branch choices for both alpha and beta to be safe.
\\
\\ Run: gp -q chlrs_forward_map_v2.gp

default(parisize, 256000000);
default(timer, 0);

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
  if(nrm == 0, error("f3inv: norm=0"));
  ni = lift(Mod(nrm, pp)^(-1));
  [((a^2 - rr*b*c)*ni) % pp, ((rr*c^2 - a*b)*ni) % pp, ((b^2 - a*c)*ni) % pp]
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

check_jac(aa, bb, t1, t2, pp) = {
  my(hh, cp, nj, target);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  target = (pp+1-t1) * (pp+1-t2);
  [nj, target, nj==target]
};

\\ cube_root_fp3: works for any pp with 3 | pp-1 and rr,rr/6 or similar being
\\ coset reps -- here we search over ALL cube roots directly instead, since
\\ we want to try all 3 branches per r explicitly (not just find one).
cube_roots_fp3(r, z3, rr, pp) = {
  my(res, c);
  res = [];
  if(r%pp==0, return([[0,0,0]]));
  if(lift(Mod(r,pp)^((pp-1)/3))==1,
    c = lift(Mod(r,pp)^( (2*(pp-1)/3+1)/3 ));  \\ placeholder, unused branch
  );
  \\ direct: find x in F_p with x^3=r (0..2 solutions if r is a cube in F_p)
  for(xx=0, pp-1, if(lift(Mod(xx,pp)^3)==r%pp, res=concat(res,[[xx,0,0]])));
  if(#res>0, return(res));
  \\ else find x with (x*a)^3 = x^3*rr = r  => x^3 = r/rr
  my(rq = lift(Mod(r,pp)*Mod(rr,pp)^(-1)));
  for(xx=0, pp-1, if(lift(Mod(xx,pp)^3)==rq, res=concat(res,[[0,xx,0]])));
  if(#res>0, return(res));
  \\ else x with (x*a^2)^3 = x^3*rr^2 = r => x^3 = r/rr^2
  my(rq2 = lift(Mod(r,pp)*Mod(rr,pp)^(-2)));
  for(xx=0, pp-1, if(lift(Mod(xx,pp)^3)==rq2, res=concat(res,[[0,0,xx]])));
  res
};

\\ ================================================================
print("=== Sanity: p=43, E1: y^2=x^3+7, E2: y^2=x^3+13 (known-good) ===");
{
  pp=43; b1=7; b2=13; rr=(-b1)%pp; z3=6;
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  t1=",t1,"  t2=",t2,"  target=",(pp+1-t1)*(pp+1-t2));
  r1=(-b1)%pp; r2=(-b2)%pp;
  avs = cube_roots_fp3(r1, z3, rr, pp);
  bvs = cube_roots_fp3(r2, z3, rr, pp);
  print("  cube roots of r1=",r1,": ",avs);
  print("  cube roots of r2=",r2,": ",bvs);
  found=0;
  for(ii=1,#avs, for(jj=1,#bvs,
    av=avs[ii]; bv=bvs[jj];
    sv=f3add(av,bv,pp); qv=f3mul(av,bv,rr,pp);
    res=richelot(sv,qv,z3,rr,pp);
    if(res[1]!=-1,
      chk=check_jac(res[1],res[2],t1,t2,pp);
      if(chk[3], print("  MATCH: alpha=",av," beta=",bv," -> a=",res[1]," b=",res[2]," #Jac=",chk[1]); found=1);
    )
  ));
  if(!found, print("  NO MATCH among ",#avs*#bvs," branch pairs"));
}
print("");

\\ ================================================================
print("=== p=1009, E1: y^2=x^3+11, E2 = quadratic twist by first non-square d ===");
{
  pp=1009; b1=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  d=2; while(kronecker(d,pp)!=-1, d++);
  b2=lift(Mod(d^3*b1,pp));
  rr=(-b1)%pp;
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  d=",d,"  b2=",b2,"  t1=",t1,"  t2=",t2,"  target=",(pp+1-t1)*(pp+1-t2));
  r1=(-b1)%pp; r2=(-b2)%pp;
  avs = cube_roots_fp3(r1, z3, rr, pp);
  bvs = cube_roots_fp3(r2, z3, rr, pp);
  print("  #cube roots of r1: ",#avs,"   #cube roots of r2: ",#bvs);
  found=0;
  for(ii=1,#avs, for(jj=1,#bvs,
    av=avs[ii]; bv=bvs[jj];
    sv=f3add(av,bv,pp); qv=f3mul(av,bv,rr,pp);
    res=richelot(sv,qv,z3,rr,pp);
    if(res[1]!=-1,
      chk=check_jac(res[1],res[2],t1,t2,pp);
      print("    branch alpha=",av," beta=",bv," -> a=",res[1]," b=",res[2]," #Jac=",chk[1]," target=",chk[2]," match=",chk[3]);
      if(chk[3], found=1);
    ,
      print("    branch alpha=",av," beta=",bv," -> cover NOT over F_p")
    )
  ));
  if(found, print("  ==> MATCH FOUND among branches"), print("  ==> NO MATCH among any of the ",#avs*#bvs," branch pairs"));
}
print("");
print("Done.");

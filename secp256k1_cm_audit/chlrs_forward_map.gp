\\ chlrs_forward_map.gp
\\
\\ Thread 3 (CHLRS Igusa forward map), autolab 2026-08-08.
\\
\\ 2026-07-27's Test 2 (howe_5pairs_v2.gp) built the naive-cover Richelot
\\ input as sv=(1+d)*alpha, qv=d*alpha^2 -- i.e. beta = d*alpha, ONE specific
\\ choice among THREE cube roots of -b2 relative to a fixed alpha (alpha is a
\\ cube root of -b1, beta must satisfy beta^3=-b2, and beta = d*z3^j*alpha for
\\ j=0,1,2 all satisfy this since z3^3=1). Only j=0 was tried and it produced
\\ #Jac=1106283 != target 1018251.
\\
\\ Hypothesis for this run: the correct Howe gluing correspondence Gamma is
\\ realized by one of the OTHER two relative twists j=1,2, not by a wholly
\\ different (undiscovered) CHLRS formula. This is a cheap, falsifiable check
\\ before committing to porting the full CHLRS Igusa-Clebsch machinery.
\\
\\ Falsifier: if none of j=0,1,2 (nor the two candidate targets #Jac in
\\ {(p+1-t1)(p+1-t2), (p+1-t1)(p+1+t2)}) match, the gluing is NOT a simple
\\ relative cube-root twist and the CHLRS Igusa-inversion route is required.
\\
\\ Run: gp -q chlrs_forward_map.gp

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
  if(nrm == 0, error("f3inv: norm=0, element not invertible"));
  ni = lift(Mod(nrm, pp)^(-1));
  [((a^2 - rr*b*c)*ni) % pp,
   ((rr*c^2 - a*b)*ni) % pp,
   ((b^2 - a*c)*ni) % pp]
};

richelot(sv, qv, z3, rr, pp) = {
  my(G1c, G1x, G2c, G2x, G3c, G3x, D0, D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,lc,lcinv,aa,bb);
  my(z3sq);
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

check_jac(aa, bb, pp) = {
  my(hh, cp, nj);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  [cp, nj]
};

\\ ================================================================
print("================================================================");
print("Sanity: p=43 baseline (j=0 reproduces howe_5pairs_v2.gp Test 1)");
print("================================================================");
{
  pp=43; rr=36; z3=6;
  sv=[0,3,0]; qv=[0,0,2];
  res=richelot(sv,qv,z3,rr,pp);
  print("  j=0: a=",res[1]," b=",res[2],"  expect a=41,b=5  -> ",
        if(res[1]==41 && res[2]==5, "MATCH", "MISMATCH"));
}
print("");

\\ ================================================================
print("================================================================");
print("Main: p=1009, E1: y^2=x^3+11, E2: y^2=x^3+d^3*11, d=first QNR.");
print("Sweep the relative cube-root twist beta = d*z3^j*alpha, j=0,1,2,");
print("against BOTH candidate targets:");
print("  T_same  = (p+1-t1)*(p+1-t2)   [Jac ~ E1 x E2]");
print("  T_twist = (p+1-t1)*(p+1+t2)   [Jac ~ E1 x E2^quadratic-twist]");
print("================================================================");
{
  pp=1009; b1=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  print("  z3 = ",z3);
  d=2; while(kronecker(d,pp)!=-1, d++);
  print("  d (first QNR) = ",d);
  b2=lift(Mod(d^3*b1,pp));
  print("  b2 = d^3*b1 mod p = ",b2);
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  t1=",t1,"  t2=",t2);
  T_same  = (pp+1-t1)*(pp+1-t2);
  T_twist = (pp+1-t1)*(pp+1+t2);
  print("  T_same=",T_same,"   T_twist=",T_twist);
  rr=(-b1)%pp;
  found=0;
  for(j=0,2,
    dj = lift(Mod(d,pp)*Mod(z3,pp)^j);
    sv=[0,(1+dj)%pp,0]; qv=[0,0,dj%pp];
    res=richelot(sv,qv,z3,rr,pp);
    if(res[1]==-1,
      print("  j=",j,": cover not defined over F_p (branch/leading-coeff failure)");
    ,
      chk=check_jac(res[1],res[2],pp);
      nj=chk[2];
      tag = if(nj==T_same, "MATCH T_same", if(nj==T_twist, "MATCH T_twist", "no match"));
      print("  j=",j,": a=",res[1]," b=",res[2],"  #Jac=",nj,"  -> ",tag);
      if(nj==T_same || nj==T_twist, found=1);
    );
  );
  if(!found,
    print("");
    print("  RESULT: none of the 3 relative cube-root twists reproduce");
    print("  Jac ~ E1 x E2 or E1 x E2^twist. The gluing correspondence is");
    print("  NOT a simple relative-cube-root choice at fixed alpha.");
  ,
    print("");
    print("  RESULT: found a working relative twist j -- forward map is");
    print("  'pick the correct cube-root branch', no new CHLRS formula needed.");
  );
}
print("");

\\ ================================================================
print("================================================================");
print("Control: does varying alpha's OWN cube-root branch (fixing beta=d*alpha)");
print("change the result, or is richelot() invariant under alpha -> z3^k*alpha?");
print("(If invariant, the j-sweep above already covers the full 3x3 = 9");
print(" combinations up to this symmetry, and the search above is complete.)");
print("================================================================");
{
  pp=1009; b1=11; d=2; while(kronecker(d,pp)!=-1, d++);
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  rr=(-b1)%pp;
  for(k=0,2,
    \\ alpha -> z3^k*alpha, beta = d*alpha (relative twist j=0 in old alpha)
    \\ sv = (1+d)*z3^k*alpha ; qv = d*z3^{2k}*alpha^2
    zk = lift(Mod(z3,pp)^k); zk2 = lift(Mod(z3,pp)^(2*k));
    sv=[0, (1+d)*zk % pp, 0]; qv=[0,0, d*zk2 % pp];
    res=richelot(sv,qv,z3,rr,pp);
    print("  alpha-branch k=",k,": a=",res[1]," b=",res[2]);
  );
}
print("");
print("Done.");

\\ ================================================================
print("================================================================");
print("Follow-up: what IS the Jacobian we got? Factor its Weil polynomial");
print("and compare its elliptic-curve factors (by trace) against E1,E2 and");
print("their quadratic twists, to identify which curve pair Richelot(naive)");
print("actually glues -- if any recognizable pair at all.");
print("================================================================");
{
  pp=1009; b1=11; d=2; while(kronecker(d,pp)!=-1, d++);
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  rr=(-b1)%pp;
  sv=[0,(1+d)%pp,0]; qv=[0,0,d%pp];
  res=richelot(sv,qv,z3,rr,pp);
  aa=res[1]; bb=res[2];
  print("  Richelot(naive cover) curve: y^2 = x^6 + ",aa,"*x^3 + ",bb," over F_",pp);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  cp = hyperellcharpoly(hh);
  print("  Weil poly P(T) = ",cp);
  fa = factor(cp);
  print("  factor(P) = ",fa);
  \\ P(T) for a genus-2 curve over F_p is T^4 - a1 T^3 + a2 T^2 - a1 p T + p^2 (real Weil poly convention varies by PARI)
  \\ Try to read it as reciprocal quartic in variable y=T (or 1/T) and recover
  \\ trace pairs (s1,s2) with s1+s2 = -coeff, s1*s2 = ... if it factors into
  \\ two quadratics over Z.
  v = variable(cp);
  print("  P(1) = #Jac = ",subst(cp,v,1));
  \\ factor() over Z already found the irreducible-factor structure above;
  \\ just read off any degree-2 factors' trace directly from their coeffs.
  found2=0;
  for(ii=1, matsize(fa)[1],
    ff = fa[ii,1];
    if(poldegree(ff)==2 && subst(ff,v,0)==pp,
      strace = -polcoeff(ff,1);
      print("  degree-2 factor T^2-(",strace,")T+",pp," (multiplicity ",fa[ii,2],")");
      found2=1;
    )
  );
  if(!found2, print("  P(T) has no degree-2 factor of the form T^2-sT+p over Z -- Jacobian is not (visibly) a product of two curves over F_p."));
  print("  For reference: t1=",pp+1-ellcard(ellinit([0,b1],pp)),
        "  t2 for E2(b2=",lift(Mod(d^3*b1,pp)),")=",pp+1-ellcard(ellinit([0,lift(Mod(d^3*b1,pp))],pp)));
}
print("Done (follow-up).");

\\ ================================================================
print("================================================================");
print("Re-check p=43 Test 1 against the CORRECT E1 x E2 target (not just E2");
print("x E2^twist, which is what howe_5pairs_v2.gp's check_jac(...,13,...)");
print("actually tested since it only takes ONE trace argument).");
print("================================================================");
{
  pp=43; b1=7; d=2; b2=lift(Mod(d^3*b1,pp));
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  E1: y^2=x^3+",b1,"  t1=",t1);
  print("  E2: y^2=x^3+",b2,"  t2=",t2,"  (b2 should be 13)");
  aa=41; bb=5; \\ the reference Richelot output
  hh=Mod(1,pp)*x^6+Mod(aa,pp)*x^3+Mod(bb,pp);
  cp=hyperellcharpoly(hh);
  print("  Weil poly = ",cp);
  print("  #Jac = ",subst(cp,variable(cp),1));
  T_same=(pp+1-t1)*(pp+1-t2);
  T_e2xtwist=(pp+1-t2)*(pp+1+t2);
  print("  T_same (E1xE2) = ",T_same,"   T_e2xtwist (E2xE2^twist) = ",T_e2xtwist);
  print("  original howe_5pairs_v2.gp check used t_expected=13 only, i.e. tested");
  print("  against T_e2xtwist, NOT T_same. Actual match target was: ",T_e2xtwist);
}
print("Done (recheck).");

\\ ================================================================
print("================================================================");
print("Is the p=43 reference Weil poly x^4+6x^3+55x^2+258x+1849 irreducible");
print("over Z too (i.e. is the p=43 'reference' case ALSO a simple Jacobian,");
print("just like the p=1009 case -- meaning the naive-cover Richelot dual has");
print("NEVER produced a split E1xE2 Jacobian, even in the baseline case)?");
print("================================================================");
{
  pp=43;
  cp43 = x^4 + 6*x^3 + 55*x^2 + 258*x + 1849;
  fa43 = factor(cp43);
  print("  factor = ",fa43);
  print("  irreducible over Z: ",matsize(fa43)[1]==1 && fa43[1,2]==1);
}
print("Done (p=43 irreducibility check).");

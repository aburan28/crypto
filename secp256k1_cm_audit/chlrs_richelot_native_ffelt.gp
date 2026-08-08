\\ chlrs_richelot_native_ffelt.gp
\\
\\ Thread 2/3 (CHLRS Igusa forward map) — Richelot-dual construction rewritten
\\ against PARI's native FFELT (ffgen) arithmetic instead of the hand-rolled
\\ F_{p^3} routines in howe_5pairs_v2.gp. Two purposes:
\\
\\ 1. Cross-validate the Z/3Z Richelot formula independently (bug check).
\\    Confirms the p=43 reference (a=41, b=5) exactly, AND is invariant under
\\    choice of alpha root / z3 root, as required by the symmetric structure
\\    of the construction (verified: relabeling alpha among its 3 conjugate
\\    roots, or z3 among its 2 roots, leaves (aa,bb) unchanged).
\\
\\ 2. Escape the old script's implicit restriction to beta = d*alpha for a
\\   *scalar* d in F_p (which only reaches pairs where b2/b1 is a perfect
\\   cube in F_p — a strict subgroup). Here beta ranges over all three
\\   genuine cube roots of -b2 in F_{p^3}, found via factor() on the FFELT
\\   polynomial ring, so any (b1,b2) pair is reachable.
\\
\\ RESULT (see chlrs_naive_cover_split_check.gp / RESEARCH_AUTOLAB_LOG.md
\\ 2026-08-08): none of the 3 beta-branches reach #Jac = #E1*#E2 for any
\\ tested pair, including the "reference" p=43 case. The Richelot-dual curve
\\ obtained this way has the SAME Frobenius zeta function (up to a full
\\ quadratic twist) as the naive cover D=(x^3+b1)(x^3+b2) itself, which
\\ chlrs_naive_cover_split_check.gp shows is never isogenous to E1 x E2.
\\ This closes out the "Z/3Z Richelot on the naive cover" sub-approach as a
\\ dead end for the forward Howe-cover construction.
\\
\\ Run: gp -q chlrs_richelot_native_ffelt.gp

default(parisize, 256000000);

richelot_native(sv, qv, z3, pp) = {
  my(z3sq, G1c,G1x,G2c,G2x,G3c,G3x, D0,D0inv);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q3r,Q6r,lc,aa,bb);
  z3sq = z3^2;
  G1c=qv;        G1x=-sv;
  G2c=z3sq*qv;   G2x=-z3*sv;
  G3c=z3*qv;     G3x=-z3sq*sv;
  D0 = (G2x*G3c-G3x*G2c) - G1x*(G3c-G2c) + G1c*(G3x-G2x);
  if (D0==0, return([-1,-1,"D0=0"]));
  D0inv = 1/D0;
  H1x2=(G3x-G2x)*D0inv; H1x1=2*(G3c-G2c)*D0inv; H1x0=(G2x*G3c-G2c*G3x)*D0inv;
  H2x2=(G1x-G3x)*D0inv; H2x1=2*(G1c-G3c)*D0inv; H2x0=(G3x*G1c-G3c*G1x)*D0inv;
  H3x2=(G2x-G1x)*D0inv; H3x1=2*(G2c-G1c)*D0inv; H3x0=(G1x*G2c-G1c*G2x)*D0inv;
  P0=H1x0*H2x0;
  P1=H1x0*H2x1+H1x1*H2x0;
  P2=H1x0*H2x2+H1x1*H2x1+H1x2*H2x0;
  P3=H1x1*H2x2+H1x2*H2x1;
  P4=H1x2*H2x2;
  Q0r=P0*H3x0;
  Q3r=P1*H3x2+P2*H3x1+P3*H3x0;
  Q6r=P4*H3x2;
  if (poldegree(Q0r.pol) > 0 || poldegree(Q3r.pol) > 0 || poldegree(Q6r.pol) > 0,
    return([-1,-1,"not-Fp"]));
  lc = Mod(polcoeff(Q6r.pol,0),pp);
  if (lc == 0, return([-1,-1,"lc=0"]));
  aa = Mod(polcoeff(Q3r.pol,0),pp)/lc;
  bb = Mod(polcoeff(Q0r.pol,0),pp)/lc;
  [lift(aa), lift(bb), "ok"];
}

check_jac(aa,bb,pp,target) = {
  my(hh,cp,nj);
  hh = Mod(1,pp)*x^6+Mod(aa,pp)*x^3+Mod(bb,pp);
  cp = hyperellcharpoly(hh);
  nj = subst(cp,variable(cp),1);
  [nj, nj==target];
}

sweep_pair(pp, b1, b2, label) = {
  my(gen,z3,rts_a,rts_b,alpha,E1,E2,t1,t2,target,found,beta,sv,qv,res,chk);
  gen = ffgen(pp^3,'a);
  z3 = Mod(lift(polrootsmod(x^2+x+1,pp)[1]),pp)*gen^0;
  rts_a = factor(x^3 - (-b1)*gen^0);
  rts_b = factor(x^3 - (-b2)*gen^0);
  alpha = -polcoeff(rts_a[1,1],0);
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  target=(pp+1-t1)*(pp+1-t2);
  print("--- ",label,": p=",pp," b1=",b1," b2=",b2,"  target(#E1*#E2)=",target," ---");
  found=0;
  for(i=1,3,
    beta = -polcoeff(rts_b[i,1],0);
    sv = alpha+beta; qv=alpha*beta;
    res = richelot_native(sv,qv,z3,pp);
    if(res[3]=="ok",
      chk=check_jac(res[1],res[2],pp,target);
      print("  beta-branch ",i,": a=",res[1]," b=",res[2],"  #Jac=",chk[1],"  MATCH=",chk[2]);
      if(chk[2],found=1);
    ,
      print("  beta-branch ",i,": ",res[3]);
    );
  );
  if(found,print("  ==> FOUND MATCHING BRANCH"),print("  ==> no branch matched #E1*#E2"));
  print("");
}

print("=== Reference check: p=43 must reproduce a=41,b=5 (scalar branch) ===");
sweep_pair(43,7,13,"p43-reference");

print("=== p=1009, genuinely non-scalar beta branches ===");
sweep_pair(1009,11,33,"p1009-b2=33");
sweep_pair(1009,11,22,"p1009-b2=22");
sweep_pair(1009,7,189,"p1009-naive(0,3)");

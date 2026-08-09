\\ chlrs_forward_search.gp
\\
\\ Thread 3 (CHLRS Igusa forward map) -- empirical search.
\\
\\ Prior finding (2026-07-27, howe_5pairs_v2.gp Test 2): using beta = d*alpha
\\ (a RATIONAL multiple of alpha, forced by taking b2 = d^3*b1) never matches
\\ #Jac(richelot dual) = #E1 * #E2 for independent E1, E2. This script tests
\\ the natural fix: use alpha, beta as INDEPENDENT cube roots of -b1, -b2 in
\\ the SAME F_{p^3} (native PARI FFELT, not hand-rolled f3* arithmetic), and
\\ sweep all root-choice combinations.
\\
\\ Richelot formula re-implemented with native FFELT ops (+,-,*,/,^) instead
\\ of the custom f3add/f3mul/f3inv library, to rule out arithmetic-transcription
\\ bugs (one was already found and fixed once in howe_5pairs.gp -> _v2.gp).

default(parisize, 256000000);

richelot_ff(sv, qv, z3) = {
  my(z3sq, G1c,G1x,G2c,G2x,G3c,G3x, D0,D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0,Q1,Q2,Q3,Q4,Q5,Q6,lc);

  z3sq = z3^2;
  G1c = qv;            G1x = -sv;
  G2c = z3sq*qv;        G2x = -(z3*sv);
  G3c = z3*qv;          G3x = -(z3sq*sv);

  D0 = (G2x*G3c - G3x*G2c) - G1x*(G3c-G2c) + G1c*(G3x-G2x);
  D0inv = D0^(-1);

  H1n2=G3x-G2x; H1n1=2*(G3c-G2c); H1n0=G2x*G3c-G2c*G3x;
  H2n2=G1x-G3x; H2n1=2*(G1c-G3c); H2n0=G3x*G1c-G3c*G1x;
  H3n2=G2x-G1x; H3n1=2*(G2c-G1c); H3n0=G1x*G2c-G1c*G2x;

  H1x2=H1n2*D0inv; H1x1=H1n1*D0inv; H1x0=H1n0*D0inv;
  H2x2=H2n2*D0inv; H2x1=H2n1*D0inv; H2x0=H2n0*D0inv;
  H3x2=H3n2*D0inv; H3x1=H3n1*D0inv; H3x0=H3n0*D0inv;

  P0=H1x0*H2x0;
  P1=H1x0*H2x1+H1x1*H2x0;
  P2=H1x0*H2x2+H1x1*H2x1+H1x2*H2x0;
  P3=H1x1*H2x2+H1x2*H2x1;
  P4=H1x2*H2x2;

  Q0=P0*H3x0;
  Q1=P0*H3x1+P1*H3x0;
  Q2=P0*H3x2+P1*H3x1+P2*H3x0;
  Q3=P1*H3x2+P2*H3x1+P3*H3x0;
  Q4=P2*H3x2+P3*H3x1+P4*H3x0;
  Q5=P3*H3x2+P4*H3x1;
  Q6=P4*H3x2;

  [Q0,Q1,Q2,Q3,Q4,Q5,Q6]
};

\\ Is an FFELT element in the base field F_p (i.e. degree-0 as a poly in the
\\ field generator)? poldegree(lift(e))<=0 for constants (or e==0).
is_base(e) = { poldegree(e.pol) <= 0 };
to_fp(e, p) = { if(e.pol==0, 0, lift(polcoef(e.pol,0))) };

\\ ================================================================
test_pair(p, b1, b2, t1, t2, label) = {
  my(T,g,z3,roots1,roots2,target,found);
  print("---- ", label, "  p=",p," b1=",b1," b2=",b2," t1=",t1," t2=",t2," ----");
  target = (p+1-t1)*(p+1-t2);
  print("  target #Jac = (p+1-t1)(p+1-t2) = ", target);

  T = ffinit(p, 3, a);
  g = ffgen(T);
  \\ primitive cube root of unity: need it in F_p if p = 1 mod 3, else in F_{p^3}
  if((p-1)%3==0,
    z3 = subst(lift(polrootsmod(x^2+x+1,p)[1]),x,1)*g^0
  ,
    z3 = polrootsmod(x^2+x+1, g)[1]
  );

  roots1 = polrootsmod(x^3 - Mod(-b1,p)*g^0, g);
  roots2 = polrootsmod(x^3 - Mod(-b2,p)*g^0, g);
  print("  #roots(-b1)=",#roots1,"  #roots(-b2)=",#roots2);

  found = 0;
  for(i=1,#roots1,
    for(j=1,#roots2,
      my(al=roots1[i], be=roots2[j], sv=al+be, qv=al*be, Q, ok, aa, bb);
      Q = richelot_ff(sv, qv, z3);
      \\ need Q1=Q2=Q4=Q5=0 identically (odd-in-x^3 vanishing) AND Q0,Q3,Q6 in F_p
      ok = (Q[2]==0 && Q[3]==0 && Q[5]==0 && Q[6]==0
            && is_base(Q[1]) && is_base(Q[4]) && is_base(Q[7]) && Q[7]!=0);
      if(ok,
        my(lc=to_fp(Q[7],p), lcinv=lift(Mod(lc,p)^(-1)));
        aa = lift(Mod(to_fp(Q[4],p)*lcinv,p));
        bb = lift(Mod(to_fp(Q[1],p)*lcinv,p));
        my(hh=Mod(1,p)*x^6+Mod(aa,p)*x^3+Mod(bb,p), cp=hyperellcharpoly(hh),
           nj=subst(cp,variable(cp),1));
        print("  (i,j)=(",i,",",j,")  a=",aa,"  b=",bb,"  #Jac=",nj,
              "  match=",nj==target);
        if(nj==target, found=1; print("    *** MATCH ***"));
      )
    )
  );
  if(!found, print("  no (i,j) combination matched target."));
  found
};

\\ ================================================================
print("================================================================");
print("Search 1: p=43 toy, E1: y^2=x^3+7 (t=?), E2: y^2=x^3+13");
print("================================================================");
{
  p=43;
  E1=ellinit([0,7],p); E2=ellinit([0,13],p);
  t1=p+1-ellcard(E1); t2=p+1-ellcard(E2);
  test_pair(p,7,13,t1,t2,"p=43 b1=7,b2=13");
}
print("");
print("================================================================");
print("Search 2: p=1009, E1: y^2=x^3+11, E2: y^2=x^3+515 (Thread2 test-2 pair)");
print("================================================================");
{
  p=1009;
  E1=ellinit([0,11],p); E2=ellinit([0,515],p);
  t1=p+1-ellcard(E1); t2=p+1-ellcard(E2);
  test_pair(p,11,515,t1,t2,"p=1009 b1=11,b2=515");
}
print("");
print("Done.");
quit

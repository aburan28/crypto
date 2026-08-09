\\ chlrs_forward_search2.gp
\\
\\ General Richelot-partition search. The 6 roots of (x^3+b1)(x^3+b2) can be
\\ split into 3 unordered pairs in 15 ways; each gives a (possibly) different
\\ Richelot-dual genus-2 curve. search1 (chlrs_forward_search.gp) only tried
\\ the 3 Galois-cyclic-compatible partitions (all giving the SAME dual up to
\\ Frobenius twist, since the summed sextic is fixed and the special zeta-
\\ twisted construction realizes only 3 of the 15 partitions). This script
\\ enumerates all 15 and checks each against the target Jacobian order.

default(parisize, 256000000);

\\ Generic Richelot dual of y^2 = G1(x)G2(x)G3(x), G_i monic quadratic
\\ x^2 + Gix*x + Gic. Returns [Q0..Q6] coeffs of the (possibly degree<6,
\\ possibly not-yet-normalized) dual polynomial.
richelot_general(G1c,G1x,G2c,G2x,G3c,G3x) = {
  my(D0,D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0,Q1,Q2,Q3,Q4,Q5,Q6);

  D0 = (G2x*G3c - G3x*G2c) - G1x*(G3c-G2c) + G1c*(G3x-G2x);
  if(D0==0, return(0));
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

is_base(e) = { poldegree(e.pol) <= 0 };
to_fp(e) = { if(e.pol==0, 0, lift(polcoef(e.pol,0))) };

\\ all 15 perfect matchings of {1,...,6} into 3 unordered pairs (hardcoded --
\\ gp's parser does not support recursive anonymous closures cleanly)
matchings6() = {
  [
   [[1,2],[3,4],[5,6]], [[1,2],[3,5],[4,6]], [[1,2],[3,6],[4,5]],
   [[1,3],[2,4],[5,6]], [[1,3],[2,5],[4,6]], [[1,3],[2,6],[4,5]],
   [[1,4],[2,3],[5,6]], [[1,4],[2,5],[3,6]], [[1,4],[2,6],[3,5]],
   [[1,5],[2,3],[4,6]], [[1,5],[2,4],[3,6]], [[1,5],[2,6],[3,4]],
   [[1,6],[2,3],[4,5]], [[1,6],[2,4],[3,5]], [[1,6],[2,5],[3,4]]
  ]
};

test_general(p, b1, b2, t1, t2, label) = {
  my(T,g,roots1,roots2,allroots,target,ms,nfound);
  print("---- ", label, "  p=",p," b1=",b1," b2=",b2," t1=",t1," t2=",t2," ----");
  target = (p+1-t1)*(p+1-t2);
  print("  target #Jac = ", target);

  T = ffinit(p,3,a); g = ffgen(T);
  roots1 = polrootsmod(x^3 - Mod(-b1,p)*g^0, g);
  roots2 = polrootsmod(x^3 - Mod(-b2,p)*g^0, g);
  allroots = concat(Vec(roots1), Vec(roots2));
  print("  6 roots gathered (3 from -b1, 3 from -b2).");

  ms = matchings6();
  print("  #matchings = ", #ms);
  nfound = 0;
  for(m=1,#ms,
    my(pr=ms[m]);
    my(r1a=allroots[pr[1][1]],r1b=allroots[pr[1][2]]);
    my(r2a=allroots[pr[2][1]],r2b=allroots[pr[2][2]]);
    my(r3a=allroots[pr[3][1]],r3b=allroots[pr[3][2]]);
    my(G1c=r1a*r1b, G1x=-(r1a+r1b));
    my(G2c=r2a*r2b, G2x=-(r2a+r2b));
    my(G3c=r3a*r3b, G3x=-(r3a+r3b));
    my(Q=richelot_general(G1c,G1x,G2c,G2x,G3c,G3x));
    if(Q==0, next());
    my(ok = Q[2]==0 && Q[3]==0 && Q[5]==0 && Q[6]==0
            && is_base(Q[1]) && is_base(Q[4]) && is_base(Q[7]) && Q[7]!=0);
    if(ok,
      my(lc=to_fp(Q[7]), lcinv=lift(Mod(lc,p)^(-1)));
      my(aa=lift(Mod(to_fp(Q[4])*lcinv,p)), bb=lift(Mod(to_fp(Q[1])*lcinv,p)));
      my(hh=Mod(1,p)*x^6+Mod(aa,p)*x^3+Mod(bb,p), cp=hyperellcharpoly(hh),
         nj=subst(cp,variable(cp),1));
      print("  matching#",m," partition=",pr,"  a=",aa," b=",bb,
            "  #Jac=",nj,"  match=",nj==target);
      if(nj==target, nfound++; print("    *** MATCH ***"));
    );
  );
  print("  total F_p-rational-descending matchings tried above; matches found = ", nfound);
  nfound
};

print("================================================================");
print("General 15-partition search, p=43 (E1: y^2=x^3+7, E2: y^2=x^3+13)");
print("================================================================");
{
  p=43; E1=ellinit([0,7],p); E2=ellinit([0,13],p);
  t1=p+1-ellcard(E1); t2=p+1-ellcard(E2);
  test_general(p,7,13,t1,t2,"p=43");
}
print("");
print("================================================================");
print("General 15-partition search, p=1009 (E1: y^2=x^3+11, E2: y^2=x^3+515)");
print("================================================================");
{
  p=1009; E1=ellinit([0,11],p); E2=ellinit([0,515],p);
  t1=p+1-ellcard(E1); t2=p+1-ellcard(E2);
  test_general(p,11,515,t1,t2,"p=1009");
}
print("Done.");
quit

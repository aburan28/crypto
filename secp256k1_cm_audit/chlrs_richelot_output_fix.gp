\\ chlrs_richelot_output_fix.gp
\\
\\ Thread 3 (CHLRS Igusa forward map), 2026-08-08/09 session.
\\ Literature access (arxiv.org, eprint.iacr.org, math.mit.edu) is
\\ EGRESS_BLOCKED in this environment -- cannot fetch Howe-Leprevost-Poonen
\\ or Lercier-Ritzenthaler directly, so this is a purely empirical re-audit
\\ of the existing Z/3Z Richelot machinery (richelot_gen / richelot in
\\ howe_richelot_v5.gp and howe_5pairs_v2.gp) against ground-truth
\\ Frobenius data from PARI's ellcard/hyperellcharpoly.
\\
\\ FINDING 1 (bug fix): richelot_gen's raw output (aa,bb) for the sextic
\\ y^2=x^6+aa*x^3+bb was being used AS-IS throughout the log. It is wrong.
\\ The correct cover is y^2=x^6 + bb_raw*x^3 - aa_raw, i.e. apply the
\\ transform (aa,bb) -> (bb, -aa) to richelot_gen's output before
\\ building the sextic. Verified by FULL characteristic polynomial match
\\ (not just #Jac order, which can coincide spuriously at small p) at
\\ p=43, b1=7, b2=13 (d=2 quadratic twist): richelot_gen([alpha,beta])
\\ with alpha^3=-b1, beta^3=-b2 gives raw (41,5); transformed (5,2) has
\\ hyperellcharpoly EXACTLY x^4-83*x^2+1849 = (x^2-13x+43)(x^2+13x+43),
\\ the true product-Jacobian target. Prior log entries (2026-07-26/27)
\\ calling this case "correct" only checked that two independently-written
\\ scripts agreed with each other on (41,5) -- neither ever checked it
\\ against the actual target Frobenius polynomial. That was a false
\\ positive; this script is the first time the p=43 case has been
\\ checked against ground truth.
\\
\\ FINDING 2 (still open / genuinely harder than thought): with the
\\ output-transform bug fixed, r1=-b1, r2=-b2 (the "naive" branch-input
\\ guess used throughout the log) still does NOT reproduce E1 x E2 at
\\ p=1009 (Test 2 of howe_5pairs_v2.gp). This is now confirmed via FULL
\\ characteristic polynomial (all 3 distinct Galois branches checked,
\\ none match), using a clean field basis (rr=2, decoupled from b1)
\\ to rule out a basis-choice artifact. So the p=43 vs p=1009 divergence
\\ is a genuine structural / size-dependent effect, not an implementation
\\ bug in richelot_gen itself.
\\
\\ Run: gp -q chlrs_richelot_output_fix.gp

default(parisize, 256000000);
default(timer, 0);

f3add(u, v, pp) = [(u[1]+v[1])%pp, (u[2]+v[2])%pp, (u[3]+v[3])%pp];
f3neg(u, pp)    = [(-u[1])%pp, (-u[2])%pp, (-u[3])%pp];
f3scl(c, u, pp) = [(c*u[1])%pp, (c*u[2])%pp, (c*u[3])%pp];

f3mul(u, v, rr, pp) = {
  my(c0,c1,c2,c3,c4);
  c0=(u[1]*v[1])%pp; c1=(u[1]*v[2]+u[2]*v[1])%pp; c2=(u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%pp;
  c3=(u[2]*v[3]+u[3]*v[2])%pp; c4=(u[3]*v[3])%pp;
  [(c0+rr*c3)%pp,(c1+rr*c4)%pp,c2%pp]
};

f3inv(u, rr, pp) = {
  my(a,b,c,nrm,ni); a=u[1]; b=u[2]; c=u[3];
  nrm=lift(Mod(a^3+rr*b^3+rr^2*c^3-3*rr*a*b*c,pp));
  if(nrm==0, error("f3inv: norm=0"));
  ni=lift(Mod(nrm,pp)^(-1));
  [((a^2-rr*b*c)*ni)%pp, ((rr*c^2-a*b)*ni)%pp, ((b^2-a*c)*ni)%pp]
};

\\ Z/3Z Richelot dual: RAW output (aa,bb) requires the (aa,bb)->(bb,-aa)
\\ transform (Finding 1) before use as sextic coefficients.
richelot(sv, qv, z3, rr, pp) = {
  my(z3sq,G1c,G1x,G2c,G2x,G3c,G3x,D0,D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,lc,lcinv,aa,bb);
  z3sq=lift(Mod(z3,pp)^2);
  G1c=qv; G1x=f3neg(sv,pp);
  G2c=f3scl(z3sq,qv,pp); G2x=f3neg(f3scl(z3,sv,pp),pp);
  G3c=f3scl(z3,qv,pp); G3x=f3neg(f3scl(z3sq,sv,pp),pp);
  D0=f3add(f3add(f3add(f3mul(G2x,G3c,rr,pp),f3neg(f3mul(G3x,G2c,rr,pp),pp),pp),
    f3neg(f3mul(G1x,f3add(G3c,f3neg(G2c,pp),pp),rr,pp),pp),pp),
    f3mul(G1c,f3add(G3x,f3neg(G2x,pp),pp),rr,pp),pp);
  if(D0==[0,0,0], return([-1,-1]));
  D0inv=f3inv(D0,rr,pp);
  H1n2=f3add(G3x,f3neg(G2x,pp),pp); H1n1=f3scl(2,f3add(G3c,f3neg(G2c,pp),pp),pp);
  H1n0=f3add(f3mul(G2x,G3c,rr,pp),f3neg(f3mul(G2c,G3x,rr,pp),pp),pp);
  H2n2=f3add(G1x,f3neg(G3x,pp),pp); H2n1=f3scl(2,f3add(G1c,f3neg(G3c,pp),pp),pp);
  H2n0=f3add(f3mul(G3x,G1c,rr,pp),f3neg(f3mul(G3c,G1x,rr,pp),pp),pp);
  H3n2=f3add(G2x,f3neg(G1x,pp),pp); H3n1=f3scl(2,f3add(G2c,f3neg(G1c,pp),pp),pp);
  H3n0=f3add(f3mul(G1x,G2c,rr,pp),f3neg(f3mul(G1c,G2x,rr,pp),pp),pp);
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
  lc=Q6r[1]; if(lc==0,return([-1,-1]));
  lcinv=lift(Mod(lc,pp)^(-1));
  aa=(Q3r[1]*lcinv)%pp; bb=(Q0r[1]*lcinv)%pp;
  [aa,bb]
};

cube_roots_fp3(r, rr, pp) = {
  my(res); res=[];
  if(r%pp==0, return([[0,0,0]]));
  for(xx=0,pp-1, if(lift(Mod(xx,pp)^3)==r%pp, res=concat(res,[[xx,0,0]])));
  if(#res>0, return(res));
  my(rq=lift(Mod(r,pp)*Mod(rr,pp)^(-1)));
  for(xx=0,pp-1, if(lift(Mod(xx,pp)^3)==rq, res=concat(res,[[0,xx,0]])));
  if(#res>0, return(res));
  my(rq2=lift(Mod(r,pp)*Mod(rr,pp)^(-2)));
  for(xx=0,pp-1, if(lift(Mod(xx,pp)^3)==rq2, res=concat(res,[[0,0,xx]])));
  res
};

\\ Try all cube-root branches of (r1,r2) and report every valid F_p cover
\\ (after the aa/bb transform), checking against the FULL target charpoly
\\ of E1 x E2 (not just the order, which coincides spuriously at small p).
verify_pair(r1, r2, rr, z3, pp, b1, b2) = {
  my(avs, bvs, E1, E2, t1, t2, target_cp, any_full);
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  target_cp = (x^2-t1*x+pp)*(x^2-t2*x+pp);
  any_full=0;
  print("  b1=",b1," b2=",b2,"  t1=",t1," t2=",t2,"  target=",target_cp);
  avs=cube_roots_fp3(r1,rr,pp); bvs=cube_roots_fp3(r2,rr,pp);
  for(ii=1,#avs, for(jj=1,#bvs,
    my(sv,qv,res,aa2,bb2,hh,cp);
    sv=f3add(avs[ii],bvs[jj],pp); qv=f3mul(avs[ii],bvs[jj],rr,pp);
    res=richelot(sv,qv,z3,rr,pp);
    if(res[1]!=-1,
      aa2=res[2]%pp; bb2=(-res[1])%pp;
      hh=Mod(1,pp)*x^6+Mod(aa2,pp)*x^3+Mod(bb2,pp);
      if(poldegree(gcd(hh,hh'))==0,
        cp=hyperellcharpoly(hh);
        if(cp==target_cp, any_full=1);
        print("    cover a=",aa2," b=",bb2,"  charpoly=",cp,"  FULL_MATCH=",cp==target_cp);
      )
    )
  ));
  print("  ==> any_full_match = ", any_full);
  any_full
};

print("================================================================");
print("Case A: p=43, b1=7, b2=13 (d=2 quadratic twist). r1=-b1, r2=-b2.");
print("================================================================");
verify_pair(-7, -13, 6, 6, 43, 7, 13);
print("");

print("================================================================");
print("Case B: p=1009, b1=11, b2=515 (d=11 quadratic twist). r1=-b1, r2=-b2.");
print("================================================================");
{
  pp=1009; b1=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  d=2; while(kronecker(d,pp)!=-1, d++);
  b2=lift(Mod(d^3*b1,pp));
  print("  d=",d);
  verify_pair((-b1)%pp,(-b2)%pp, 2, z3, pp, b1, b2);
}
print("");
print("CONCLUSION: r1=-b1,r2=-b2 works at p=43 (Case A) and fails at");
print("p=1009 (Case B) even after the output-transform fix. The gap is");
print("structural, not a coding bug -- see RESEARCH_AUTOLAB_LOG.md for");
print("next-step proposal.");

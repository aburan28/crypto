\\ chlrs_forward_map.gp
\\
\\ Thread 3 (CHLRS Igusa forward map), proposed 2026-07-27.
\\ Reuses the Z/3Z Richelot machinery from howe_5pairs_v2.gp (verified correct
\\ on p=43) and tests the hypothesis that the 2026-07-27 "Test 2" failure
\\ (p=1009, pair with d=first-non-square) was caused by picking the WRONG
\\ F_{p^3}-conjugate of beta, i.e. the wrong Galois-equivariant bijection
\\ E1[2] -> E2[2], rather than by a missing "Igusa inversion" formula.
\\
\\ E1: y^2=x^3+b1, 2-torsion x-coords are the 3 roots of x^3=-b1: alpha,
\\ z3*alpha, z3^2*alpha (alpha in F_{p^3} when -b1 is not a cube in F_p).
\\ E2: y^2=x^3+b2, roots beta, z3*beta, z3^2*beta similarly.
\\
\\ A Galois-equivariant bijection between the two 2-torsion orbits must send
\\ alpha -> beta_k := z3^k * beta for SOME FIXED k in {0,1,2} (all 3 conjugates
\\ of alpha map correspondingly, since Frobenius commutes with mult-by-z3^k).
\\ The existing script only tried k=0. This script tries k=0,1,2 and checks
\\ which (if any) reproduces #Jac(C) = #E1 * #E2.
\\
\\ Run: gp -q chlrs_forward_map.gp

default(parisize, 256000000);
default(timer, 0);

\\ ================================================================
\\ F_{pp^3} arithmetic (verbatim from howe_5pairs_v2.gp)
\\ ================================================================

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

\\ Given alpha (F_{p^3} elt, root of x^3=rr1=-b1) and a candidate beta
\\ (F_{p^3} elt, root of x^3=rr1... no, root of x^3=-b2), run richelot and
\\ report whether #Jac(C) matches #E1*#E2.
try_pairing(alpha, beta, z3, rr1, pp, t1, t2, label) = {
  my(sv, qv, res, chk);
  sv = f3add(alpha, beta, pp);
  qv = f3mul(alpha, beta, rr1, pp);
  res = richelot(sv, qv, z3, rr1, pp);
  if(res[1]==-1,
    print("  [",label,"] cover NOT over F_p (Richelot denominator/leading coeff issue)");
    return(0)
  );
  chk = check_jac(res[1], res[2], t1, pp);
  print("  [",label,"] a=",res[1],"  b=",res[2],"  #Jac=",chk[2],
        "  target(#E1*#E2)=",chk[3],"  MATCH=",chk[4]);
  chk[4]
};

\\ ================================================================
print("================================================================");
print("Regression: p=43 pair from howe_5pairs_v2.gp Test 1 (k=0 known-good)");
print("================================================================");
{
  pp=43; rr=36; z3=6; z3sq=lift(Mod(z3,pp)^2);
  alpha=[0,1,0]; \\ alpha itself (b1=7 => rr=-7=36 mod 43, alpha^3=rr)
  beta_k0=[0,3,0]; \\ matches sv=[0,3,0] when alpha=[0,1,0] i.e. beta=[0,2,0]... reuse original test directly
  \\ Original test used sv=[0,3,0], qv=[0,0,2] directly (not alpha,beta explicitly).
  \\ Reconstruct: alpha=[0,1,0] (coefficient 1), beta s.t. alpha+beta=[0,3,0] => beta=[0,2,0].
  \\ Check alpha*beta should be qv=[0,0,2]: f3mul([0,1,0],[0,2,0],36,43) should be [0,0,2].
  print("  alpha*beta check: ", f3mul([0,1,0],[0,2,0],36,43), " (expect [0,0,2])");
  res=try_pairing([0,1,0],[0,2,0],z3,rr,pp,13,13,"k=0 (regression)");
}
print("");

\\ ================================================================
print("================================================================");
print("Thread 3: p=1009, E1: b1=11, E2: b2=d^3*b1 for first non-square d.");
print("Testing all 3 Galois-conjugate pairings alpha <-> z3^k * beta.");
print("================================================================");
{
  pp=1009; b1=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]); z3sq=lift(Mod(z3,pp)^2);
  d=2; while(kronecker(d,pp)!=-1, d++);
  b2=lift(Mod(d^3*b1,pp));
  print("  z3=",z3,"  d=",d,"  b2=",b2);
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  trace E1=",t1,"  trace E2=",t2,"  #E1*#E2=",ellcard(E1)*ellcard(E2));
  rr1=(-b1)%pp;
  alpha=[0,1,0]; \\ alpha = the F_{p^3} generator itself: alpha^3 = rr1 = -b1
  beta0=[0,d%pp,0];                 \\ beta = d*alpha        (k=0, the 07-27 attempt)
  beta1=[0,lift(Mod(d*z3,pp)),0];   \\ beta = d*z3*alpha     (k=1)
  beta2=[0,lift(Mod(d*z3sq,pp)),0]; \\ beta = d*z3^2*alpha   (k=2)
  any_match = 0;
  any_match += try_pairing(alpha, beta0, z3, rr1, pp, t1, t2, "k=0");
  any_match += try_pairing(alpha, beta1, z3, rr1, pp, t1, t2, "k=1");
  any_match += try_pairing(alpha, beta2, z3, rr1, pp, t1, t2, "k=2");
  print("");
  if(any_match>0,
    print("  RESULT: forward map found -- Howe cover = Richelot(alpha, z3^k*beta) for some k."),
    print("  RESULT: none of the 3 conjugate pairings reproduce #E1*#E2.");
    print("  Also trying #E1 * #E2^twist and #E1^twist * #E2 as alternate targets:");
    for(kk=0,2,
      bb = if(kk==0,beta0,if(kk==1,beta1,beta2));
      sv=f3add(alpha,bb,pp); qv=f3mul(alpha,bb,rr1,pp);
      res=richelot(sv,qv,z3,rr1,pp);
      if(res[1]!=-1,
        hh = Mod(1,pp)*x^6 + Mod(res[1],pp)*x^3 + Mod(res[2],pp);
        cp = hyperellcharpoly(hh); nj = subst(cp,variable(cp),1);
        print("    k=",kk,": #Jac=",nj,
              "  E1*E2tw=",(pp+1-t1)*(pp+1+t2),
              "  E1tw*E2=",(pp+1+t1)*(pp+1-t2),
              "  E1tw*E2tw=",(pp+1+t1)*(pp+1+t2))
      )
    )
  );
}
print("");
print("================================================================");
print("Done.");
print("================================================================");

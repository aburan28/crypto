\\ chlrs_bb_invariant_general.gp
\\
\\ Generalises chlrs_bb_invariant.gp (which hardcoded p=43, rr=36) to
\\ arbitrary p (p=1 mod 3) and arbitrary non-cube rr, to check the
\\ bb_dual = r1*r2 mod p invariant is not a p=43 coincidence.
\\ GP 2.15.4 does not support nested {} blocks inside function bodies, so
\\ p/rr/z3 are globals re-set before each run_test() call (same style as
\\ the other scripts in this directory).
\\
\\ Run: gp -q chlrs_bb_invariant_general.gp

default(parisize, 256000000);
default(timer, 0);

fadd(u,v) = [(u[1]+v[1])%p,(u[2]+v[2])%p,(u[3]+v[3])%p];
fneg(u) = [(-u[1])%p,(-u[2])%p,(-u[3])%p];
fscl(c,u) = [(c*u[1])%p,(c*u[2])%p,(c*u[3])%p];
fmul(u,v) = [(u[1]*v[1]+rr*(u[2]*v[3]+u[3]*v[2]))%p,(u[1]*v[2]+u[2]*v[1]+rr*u[3]*v[3])%p,(u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p];
finv(u) = {my(a,b,c,n,ni); a=u[1]; b=u[2]; c=u[3]; n=lift(Mod(a^3+rr*b^3+rr^2*c^3-3*rr*a*b*c,p)); if(n==0,error("f")); ni=lift(Mod(n,p)^(-1)); [(a^2-rr*b*c)*ni%p,(rr*c^2-a*b)*ni%p,(b^2-a*c)*ni%p]};
is_cube(r) = (lift(Mod(r,p)^((p-1)/3))==1);
cube_root_fp3(r) = {my(r1,r2,r3,c); r=r%p; if(r==0,return([0,0,0]));
  r1=r; if(is_cube(r1),c=lift(Mod(r1,p)^((2*p-1)/3));return([c,0,0]));
  r2=(r*lift(Mod(rr,p)^(-1)))%p; if(is_cube(r2),c=lift(Mod(r2,p)^((2*p-1)/3));return([0,c,0]));
  r3=(r*lift(Mod(rr,p)^(-2)))%p; if(is_cube(r3),c=lift(Mod(r3,p)^((2*p-1)/3));return([0,0,c]));
  error("cube_root_fp3 failed for r=",r," p=",p," rr=",rr)};

richelot_gen(sv,qv) = {
  my(G1c,G1x,G2c,G2x,G3c,G3x,D0,D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,lc,lcinv,aa,bb);
  G1c=qv; G1x=fneg(sv);
  G2c=fscl(z3sq,qv); G2x=fneg(fscl(z3,sv));
  G3c=fscl(z3,qv);   G3x=fneg(fscl(z3sq,sv));
  D0=fadd(fadd(fadd(fmul(G2x,G3c),fneg(fmul(G3x,G2c))),fneg(fmul(G1x,fadd(G3c,fneg(G2c))))),fmul(G1c,fadd(G3x,fneg(G2x))));
  D0inv=finv(D0);
  H1n2=fadd(G3x,fneg(G2x)); H1n1=fscl(2,fadd(G3c,fneg(G2c))); H1n0=fadd(fmul(G2x,G3c),fneg(fmul(G2c,G3x)));
  H2n2=fadd(G1x,fneg(G3x)); H2n1=fscl(2,fadd(G1c,fneg(G3c))); H2n0=fadd(fmul(G3x,G1c),fneg(fmul(G3c,G1x)));
  H3n2=fadd(G2x,fneg(G1x)); H3n1=fscl(2,fadd(G2c,fneg(G1c))); H3n0=fadd(fmul(G1x,G2c),fneg(fmul(G1c,G2x)));
  H1x2=fmul(H1n2,D0inv); H1x1=fmul(H1n1,D0inv); H1x0=fmul(H1n0,D0inv);
  H2x2=fmul(H2n2,D0inv); H2x1=fmul(H2n1,D0inv); H2x0=fmul(H2n0,D0inv);
  H3x2=fmul(H3n2,D0inv); H3x1=fmul(H3n1,D0inv); H3x0=fmul(H3n0,D0inv);
  P0=fmul(H1x0,H2x0); P1=fadd(fmul(H1x0,H2x1),fmul(H1x1,H2x0));
  P2=fadd(fadd(fmul(H1x0,H2x2),fmul(H1x1,H2x1)),fmul(H1x2,H2x0));
  P3=fadd(fmul(H1x1,H2x2),fmul(H1x2,H2x1)); P4=fmul(H1x2,H2x2);
  Q0r=fmul(P0,H3x0); Q1r=fadd(fmul(P0,H3x1),fmul(P1,H3x0));
  Q2r=fadd(fadd(fmul(P0,H3x2),fmul(P1,H3x1)),fmul(P2,H3x0));
  Q3r=fadd(fadd(fmul(P1,H3x2),fmul(P2,H3x1)),fmul(P3,H3x0));
  Q4r=fadd(fadd(fmul(P2,H3x2),fmul(P3,H3x1)),fmul(P4,H3x0));
  Q5r=fadd(fmul(P3,H3x2),fmul(P4,H3x1)); Q6r=fmul(P4,H3x2);
  if(Q0r[2]!=0||Q0r[3]!=0||Q3r[2]!=0||Q3r[3]!=0||Q6r[2]!=0||Q6r[3]!=0,return([-1,-1]));
  lc=Q6r[1]; if(lc==0,return([-1,-1]));
  lcinv=lift(Mod(lc,p)^(-1));
  aa=(Q3r[1]*lcinv)%p; bb=(Q0r[1]*lcinv)%p;
  [aa,bb]};

run_test(pp) = {
  p = pp;
  z3 = lift(polrootsmod(x^2+x+1,p)[1]); z3sq = lift(Mod(z3,p)^2);
  rr = 2; while(lift(Mod(rr,p)^((p-1)/3))==1, rr++);
  ntest=0; nok=0; nmismatch=0;
  for(r1=1,p-1,
    av=cube_root_fp3(r1);
    for(r2=1,p-1,
      bv=cube_root_fp3(r2);
      for(kk=0,2,
        bvk=if(kk==0,bv,if(kk==1,fscl(z3,bv),fscl(z3sq,bv)));
        s=fadd(av,bvk); q=fmul(av,bvk);
        res=trap(,[-1,-1],richelot_gen(s,q));
        if(res[1]!=-1,
          ntest++;
          expect=(r1*r2)%p;
          if(res[2]==expect, nok++, nmismatch++)
        )
      )
    )
  );
  print("  p=",p,"  rr=",rr,"  bb==r1*r2:  ",nok,"/",ntest," match  (",nmismatch," mismatches)");
};

print("Cross-prime check of bb_dual = r1*r2 invariant:");
run_test(43);
run_test(61);
run_test(67);
run_test(79);
print("Done.");

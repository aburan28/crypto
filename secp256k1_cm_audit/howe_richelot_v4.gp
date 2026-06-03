\\ howe_richelot_v4.gp
\\ For each of the 7 Z/3Z canonical classes over F_43, compute the
\\ Richelot image and check which class maps to the target isogeny class.
\\
\\ Key insight: the Howe-glued curve C satisfies: the Richelot of Jac(C)
\\ under sigma_0 gives another curve ALSO in the target isogeny class
\\ (T^4-83T^2+1849). By tracking the Z/3Z Richelot graph we identify C.
\\
\\ Additionally: compute Igusa invariants of all 7 classes via transvectants
\\ to get a fingerprint.
\\
\\ Run: gp -q howe_richelot_v4.gp

default(parisize, 64000000); default(timer, 0);
p = 43; z3 = 6; z3sq = 36;

\\ F_{43^3} = F_43[t]/(t^3 - 36). Basis: [1, t, t^2]. t^3=36.
fadd(u,v) = [(u[1]+v[1])%p,(u[2]+v[2])%p,(u[3]+v[3])%p];
fneg(u) = [(-u[1])%p,(-u[2])%p,(-u[3])%p];
fscl(c,u) = [(c*u[1])%p,(c*u[2])%p,(c*u[3])%p];
fcon(c) = [c%p,0,0];
fmul(u,v) = [(u[1]*v[1]+36*(u[2]*v[3]+u[3]*v[2]))%p,(u[1]*v[2]+u[2]*v[1]+36*u[3]*v[3])%p,(u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p];
finv(u) = {my(a,b,c,n); a=u[1]; b=u[2]; c=u[3];
  n = lift(Mod(a^3+36*b^3+36^2*c^3-3*36*a*b*c, p));
  if(n==0, error("non-invertible")); ni=lift(Mod(n,p)^(-1));
  [(a^2-36*b*c)*ni%p, (36*c^2-a*b)*ni%p, (b^2-a*c)*ni%p]};

\\ Cube roots in F_{43^3}:
\\ Any r in F_43* can be written as c^3, c^3*36, or c^3*6.
\\ The cube root "c" satisfying c^3 = r:
\\ If r is a cube (r^14≡1 mod 43): r^{1/3} = r^5 mod 43.
\\ Algorithm: find alpha in F_{43^3} with alpha^3=r by trying:
\\   alpha = c (scalar): c^3=r, if r a cube
\\   alpha = c*t:       c^3*36=r, i.e., c^3=r/36, if r/36 a cube
\\   alpha = c*t^2:     c^3*6=r,  i.e., c^3=r/6, if r/6 a cube
is_cube(r) = (lift(Mod(r,p)^((p-1)/3)) == 1);
cube_root(r) = lift(Mod(r,p)^5);  \\ r^5 = r^{1/3} if r is a cube

cube_root_fp3(r) = {
  my(r1,r2,r3,c);
  r = r%p; if(r==0, return([0,0,0]));
  \\ Try alpha = c (scalar, c^3=r)
  r1 = r;
  if(is_cube(r1), c=cube_root(r1); return([c,0,0]));
  \\ Try alpha = c*t (c^3=r/36)
  r2 = (r*lift(Mod(36,p)^(-1)))%p;
  if(is_cube(r2), c=cube_root(r2); return([0,c,0]));
  \\ Try alpha = c*t^2 (c^3=r/6)
  r3 = (r*lift(Mod(6,p)^(-1)))%p;
  if(is_cube(r3), c=cube_root(r3); return([0,0,c]));
  error("cube_root_fp3 failed for r=",r)
};

\\ Verify cube_root_fp3
test = cube_root_fp3(36); tt3 = fmul(fmul(test,test),test);
print("cube_root_fp3(36) = ",test,"  cubed = ",tt3," (expect [36,0,0])");
test2 = cube_root_fp3(30); tt3_2 = fmul(fmul(test2,test2),test2);
print("cube_root_fp3(30) = ",test2,"  cubed = ",tt3_2," (expect [30,0,0])");
test3 = cube_root_fp3(5);  tt3_3 = fmul(fmul(test3,test3),test3);
print("cube_root_fp3(5)  = ",test3,"  cubed = ",tt3_3," (expect [5,0,0])");
test4 = cube_root_fp3(24); tt3_4 = fmul(fmul(test4,test4),test4);
print("cube_root_fp3(24) = ",test4,"  cubed = ",tt3_4," (expect [24,0,0])");
print("");

\\ General Z/3Z Richelot with s,q as F_{43^3} vectors.
\\ Returns monic [a,b] of image y^2=x^6+a*x^3+b, or [-1,-1] on error.
richelot_gen(sv, qv) = {
  my(G1c,G1x,G2c,G2x,G3c,G3x,Gx2,D0,D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,lc,lcinv,aa,bb);
  Gx2=fcon(1);
  G1c=qv; G1x=fneg(sv);
  G2c=fscl(z3sq,qv); G2x=fneg(fscl(z3,sv));
  G3c=fscl(z3,qv); G3x=fneg(fscl(z3sq,sv));
  \\ Delta = det 3x3
  D0=fadd(fadd(fmul(Gx2,fadd(fmul(G2x,G3c),fneg(fmul(G3x,G2c)))),fneg(fmul(G1x,fadd(fmul(Gx2,G3c),fneg(fmul(Gx2,G2c)))))),fmul(G1c,fadd(fmul(Gx2,G3x),fneg(fmul(Gx2,G2x)))));
  \\ H1 = (G2'G3 - G2G3') / D0,  G_i' = 2x + G_ix
  H1n2=fadd(fadd(fmul(fcon(2),G3x),fmul(G2x,Gx2)),fneg(fadd(fmul(Gx2,G3x),fmul(G2x,fcon(2)))));
  H1n1=fadd(fadd(fmul(fcon(2),G3c),fmul(G2x,Gx2)),fneg(fadd(fmul(G2x,G3x),fmul(G2c,fcon(2)))));
  H1n0=fadd(fmul(G2x,G3c),fneg(fmul(G2c,G3x)));
  H2n2=fadd(fadd(fmul(fcon(2),G1x),fmul(G3x,Gx2)),fneg(fadd(fmul(Gx2,G1x),fmul(G3x,fcon(2)))));
  H2n1=fadd(fadd(fmul(fcon(2),G1c),fmul(G3x,Gx2)),fneg(fadd(fmul(G3x,G1x),fmul(G3c,fcon(2)))));
  H2n0=fadd(fmul(G3x,G1c),fneg(fmul(G3c,G1x)));
  H3n2=fadd(fadd(fmul(fcon(2),G2x),fmul(G1x,Gx2)),fneg(fadd(fmul(Gx2,G2x),fmul(G1x,fcon(2)))));
  H3n1=fadd(fadd(fmul(fcon(2),G2c),fmul(G1x,Gx2)),fneg(fadd(fmul(G1x,G2x),fmul(G1c,fcon(2)))));
  H3n0=fadd(fmul(G1x,G2c),fneg(fmul(G1c,G2x)));
  D0inv=finv(D0);
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
  \\ Check all in F_43
  if(Q0r[2]!=0||Q0r[3]!=0||Q3r[2]!=0||Q3r[3]!=0||Q6r[2]!=0||Q6r[3]!=0, return([-1,-1]));
  lc=Q6r[1]; if(lc==0, return([-1,-1]));
  lcinv=lift(Mod(lc,p)^(-1));
  aa=(Q3r[1]*lcinv)%p; bb=(Q0r[1]*lcinv)%p;
  [aa,bb]
};

\\ Verify richelot_gen matches richelot_zt3 for naive cover (S=3,Q=2)
sv0=fscl(3,[0,1,0]); qv0=fscl(2,[0,0,1]);
res0=richelot_gen(sv0,qv0);
print("Naive cover sigma_0 check: a=",res0[1],"  b=",res0[2],"  (expect a=41 b=5)");
print("");

\\ ================================================================
\\ For each of the 7 canonical Z/3Z classes, find cube roots and
\\ compute Richelot images for all 3 matchings (sigma_0, sigma_1, sigma_2).
\\ ================================================================
print("================================================================");
print("Z/3Z Richelot graph for 7 canonical classes over F_43");
print("================================================================");
print("");

ca=[1,2,3,4,6,8,16]; cb=[2,8,42,32,39,42,39];
two_inv = lift(Mod(2,p)^(-1));

for (ci=1,7,
  aa=ca[ci]; bb=cb[ci];
  print("--- Class ",ci,": a=",aa,"  b=",bb,"  (y^2=x^6+",aa,"x^3+",bb,") ---");
  disc_ab = (aa^2 - 4*bb) % p;
  \\ Check if disc is a square
  sq_test = lift(Mod(disc_ab,p)^((p-1)/2));
  if (sq_test == 1,
    sqrtd = lift(Mod(disc_ab,p)^((p+1)/4));
    r1 = ((-aa + sqrtd + 4*p) * two_inv) % p;
    r2 = ((-aa - sqrtd + 4*p) * two_inv) % p;
    print("  disc=",disc_ab,"  sqrt=",sqrtd,"  r1=",r1,"  r2=",r2);
    \\ Find cube roots
    alpha_v = cube_root_fp3(r1);
    beta_v = cube_root_fp3(r2);
    print("  alpha^3=r1: alpha=",alpha_v,"  beta^3=r2: beta=",beta_v);
    \\ Verify
    a3 = fmul(fmul(alpha_v,alpha_v),alpha_v);
    b3 = fmul(fmul(beta_v,beta_v),beta_v);
    print("  alpha^3=",a3,"  beta^3=",b3,"  (r1=[",r1,",0,0] r2=[",r2,",0,0])");
    \\ Compute sigma_0: pair (alpha,beta),(z3*alpha,z3*beta),(z3^2*alpha,z3^2*beta)
    s0v = fadd(alpha_v, beta_v);
    q0v = fmul(alpha_v, beta_v);
    res_s0 = richelot_gen(s0v, q0v);
    \\ sigma_1: pair (alpha,z3*beta),(z3*alpha,z3^2*beta),(z3^2*alpha,beta)
    s1v = fadd(alpha_v, fscl(z3,beta_v));
    q1v = fmul(alpha_v, fscl(z3,beta_v));
    res_s1 = richelot_gen(s1v, q1v);
    \\ sigma_2: pair (alpha,z3^2*beta),(z3*alpha,beta),(z3^2*alpha,z3*beta)
    s2v = fadd(alpha_v, fscl(z3sq,beta_v));
    q2v = fmul(alpha_v, fscl(z3sq,beta_v));
    res_s2 = richelot_gen(s2v, q2v);
    print("  Richelot sigma_0: a=",res_s0[1],"  b=",res_s0[2]);
    print("  Richelot sigma_1: a=",res_s1[1],"  b=",res_s1[2]);
    print("  Richelot sigma_2: a=",res_s2[1],"  b=",res_s2[2]);
    \\ Verify char polys
    for (si=0,2,
      rv=[res_s0,res_s1,res_s2][si+1];
      if(rv[1]>=0,
        h=Mod(1,p)*x^6+Mod(rv[1],p)*x^3+Mod(rv[2],p);
        cp=hyperellcharpoly(h); nj=subst(cp,variable(cp),1);
        print("    sigma_",si," char poly: ",cp,"  #Jac=",nj)
      )
    )
  ,
    if (disc_ab == 0,
      print("  disc=0 (r1=r2): degenerate")
    ,
      print("  disc=",disc_ab," not a square in F_43 => r1,r2 in F_{43^2}")
      \\ Skip for now: handling F_{43^2} roots requires additional work
    )
  );
  print("")
);

print("================================================================");
print("Note: 'MATCH' expected when image has char poly T^4-83T^2+1849 and");
print("  a=target_class_a, b=target_class_b.");
print("Howe-glued class = the one whose sigma_k image is ALSO in {7 classes}");
print("AND whose img connects to the E1 x E2 structure.");
print("================================================================");

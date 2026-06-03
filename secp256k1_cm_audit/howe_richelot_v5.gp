\\ howe_richelot_v5.gp
\\ Fix two bugs from v4:
\\   1. H_in1 formula used Gx2(=1) instead of G_jx*G_kx — fixed
\\   2. Top-level for-loop with multi-line body — moved into helper fn do_class()
\\ Run: gp -q howe_richelot_v5.gp

default(parisize, 64000000); default(timer, 0);
p = 43; z3 = 6; z3sq = 36;

fadd(u,v) = [(u[1]+v[1])%p,(u[2]+v[2])%p,(u[3]+v[3])%p];
fneg(u) = [(-u[1])%p,(-u[2])%p,(-u[3])%p];
fscl(c,u) = [(c*u[1])%p,(c*u[2])%p,(c*u[3])%p];
fcon(c) = [c%p,0,0];
fmul(u,v) = [(u[1]*v[1]+36*(u[2]*v[3]+u[3]*v[2]))%p,(u[1]*v[2]+u[2]*v[1]+36*u[3]*v[3])%p,(u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p];
finv(u) = {my(a,b,c,n,ni); a=u[1]; b=u[2]; c=u[3]; n=lift(Mod(a^3+36*b^3+36^2*c^3-3*36*a*b*c,p)); if(n==0,error("finv: norm=0")); ni=lift(Mod(n,p)^(-1)); [(a^2-36*b*c)*ni%p,(36*c^2-a*b)*ni%p,(b^2-a*c)*ni%p]};

is_cube(r) = (lift(Mod(r,p)^((p-1)/3))==1);
cube_root_fp3(r) = {my(r1,r2,r3,c); r=r%p; if(r==0,return([0,0,0])); r1=r; if(is_cube(r1),c=lift(Mod(r1,p)^5);return([c,0,0])); r2=(r*lift(Mod(36,p)^(-1)))%p; if(is_cube(r2),c=lift(Mod(r2,p)^5);return([0,c,0])); r3=(r*lift(Mod(6,p)^(-1)))%p; if(is_cube(r3),c=lift(Mod(r3,p)^5);return([0,0,c])); error("cube_root_fp3 failed for r=",r)};

\\ General Z/3Z Richelot: G_i = x^2 + G_ix*x + G_ic
\\ H_i = (G_j'*G_k - G_j*G_k') / Delta  where G'(x)=2x+G_ix
\\ H_in2 = G_kx - G_jx
\\ H_in1 = 2*(G_kc - G_jc)          [FIXED from v4: was G_jx*1, now 2*(Gkc-Gjc)]
\\ H_in0 = G_jx*G_kc - G_jc*G_kx
richelot_gen(sv,qv) = {
  my(G1c,G1x,G2c,G2x,G3c,G3x,D0,D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,lc,lcinv,aa,bb);
  G1c=qv; G1x=fneg(sv);
  G2c=fscl(z3sq,qv); G2x=fneg(fscl(z3,sv));
  G3c=fscl(z3,qv);   G3x=fneg(fscl(z3sq,sv));
  \\ Delta = G1x*(G2x*G3c - G3x*G2c) - but using 3x3 det with Gx2=1:
  \\ Delta = 1*(G2x*G3c - G3x*G2c) - G1x*(1*G3c - 1*G2c) + G1c*(1*G3x - 1*G2x)
  D0=fadd(fadd(fadd(fmul(G2x,G3c),fneg(fmul(G3x,G2c))),fneg(fmul(G1x,fadd(G3c,fneg(G2c))))),fmul(G1c,fadd(G3x,fneg(G2x))));
  D0inv=finv(D0);
  \\ H1 numerator (cyclic: j=2,k=3)
  H1n2=fadd(G3x,fneg(G2x));
  H1n1=fscl(2,fadd(G3c,fneg(G2c)));
  H1n0=fadd(fmul(G2x,G3c),fneg(fmul(G2c,G3x)));
  \\ H2 numerator (cyclic: j=3,k=1)
  H2n2=fadd(G1x,fneg(G3x));
  H2n1=fscl(2,fadd(G1c,fneg(G3c)));
  H2n0=fadd(fmul(G3x,G1c),fneg(fmul(G3c,G1x)));
  \\ H3 numerator (cyclic: j=1,k=2)
  H3n2=fadd(G2x,fneg(G1x));
  H3n1=fscl(2,fadd(G2c,fneg(G1c)));
  H3n0=fadd(fmul(G1x,G2c),fneg(fmul(G1c,G2x)));
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

\\ Verify sigma_0 of naive cover: sv=[0,3,0], qv=[0,0,2] => expect a=41, b=5
res_check = richelot_gen([0,3,0],[0,0,2]);
print("Naive cover sigma_0: a=",res_check[1],"  b=",res_check[2],"  (expect a=41 b=5)");

\\ Helper: compute charply #Jac for y^2=x^6+aa*x^3+bb
charpoly_jac(aa,bb) = {my(h,cp,nj); h=Mod(1,p)*x^6+Mod(aa,p)*x^3+Mod(bb,p); cp=hyperellcharpoly(h); nj=subst(cp,variable(cp),1); [cp,nj]};

\\ Helper: process one canonical class
do_class(ci,ca,cb) = {
  my(aa,bb,disc_ab,sq_test,sqrtd,two_inv,r1,r2,av,bv,a3,b3);
  my(s0v,q0v,s1v,q1v,s2v,q2v,res0,res1,res2,cpd,tag);
  aa=ca[ci]; bb=cb[ci];
  print("--- Class ",ci,": a=",aa,"  b=",bb," ---");
  disc_ab=(aa^2-4*bb)%p;
  sq_test=lift(Mod(disc_ab,p)^((p-1)/2));
  if(sq_test!=1, print("  disc=",disc_ab," not a square => skip"); print(""); return());
  sqrtd=lift(Mod(disc_ab,p)^((p+1)/4));
  two_inv=lift(Mod(2,p)^(-1));
  r1=((-aa+sqrtd+4*p)*two_inv)%p;
  r2=((-aa-sqrtd+4*p)*two_inv)%p;
  print("  disc=",disc_ab," sqrtd=",sqrtd," r1=",r1," r2=",r2);
  av=cube_root_fp3(r1); bv=cube_root_fp3(r2);
  a3=fmul(fmul(av,av),av); b3=fmul(fmul(bv,bv),bv);
  print("  alpha=",av," alpha^3=",a3,"  beta=",bv," beta^3=",b3);
  s0v=fadd(av,bv);                     q0v=fmul(av,bv);
  s1v=fadd(av,fscl(z3,bv));            q1v=fmul(av,fscl(z3,bv));
  s2v=fadd(av,fscl(z3sq,bv));          q2v=fmul(av,fscl(z3sq,bv));
  res0=richelot_gen(s0v,q0v);
  res1=richelot_gen(s1v,q1v);
  res2=richelot_gen(s2v,q2v);
  print("  sigma_0: a=",res0[1],"  b=",res0[2]);
  print("  sigma_1: a=",res1[1],"  b=",res1[2]);
  print("  sigma_2: a=",res2[1],"  b=",res2[2]);
  if(res0[1]>=0, cpd=charpoly_jac(res0[1],res0[2]); print("    sigma_0 cp=",cpd[1],"  #Jac=",cpd[2]));
  if(res1[1]>=0, cpd=charpoly_jac(res1[1],res1[2]); print("    sigma_1 cp=",cpd[1],"  #Jac=",cpd[2]));
  if(res2[1]>=0, cpd=charpoly_jac(res2[1],res2[2]); print("    sigma_2 cp=",cpd[1],"  #Jac=",cpd[2]));
  print("")};

ca=[1,2,3,4,6,8,16]; cb=[2,8,42,32,39,42,39];
print("");
print("================================================================");
print("Z/3Z Richelot graph for 7 canonical classes over F_43");
print("Target char poly: T^4-83*T^2+1849 = (T^2-13T+43)(T^2+13T+43)");
print("================================================================");
print("");
do_class(1,ca,cb);
do_class(2,ca,cb);
do_class(3,ca,cb);
do_class(4,ca,cb);
do_class(5,ca,cb);
do_class(6,ca,cb);
do_class(7,ca,cb);
print("================================================================");
print("Done.");

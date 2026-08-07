\\ chlrs_bruteforce_charpoly.gp
\\
\\ Strengthens chlrs_bruteforce_r1r2.gp: that script's "hits" only matched
\\ #Jac(F_p) (a single integer), which is a WEAK test -- genus-2 Jacobian
\\ orders over F_p=43 live in a narrow Hasse-Weil window (~1700 valid covers
\\ found, likely far fewer *distinct* orders), so matching the scalar count
\\ is expected to produce false positives by pure collision.
\\
\\ This script requires the FULL degree-4 Frobenius characteristic
\\ polynomial of Jac(C) to equal (T^2-t1*T+p)*(T^2-t2*T+p) exactly -- the
\\ real Honda-Tate statement of "Jac(C) ~ E1 x E2" -- not just matching
\\ order. Also reports how many DISTINCT char polys occur among all valid
\\ covers, to quantify how weak the scalar-count test was.
\\
\\ Run: gp -q chlrs_bruteforce_charpoly.gp

default(parisize, 256000000);
default(timer, 0);

p = 43; z3 = 6; z3sq = 36; rr = 36;

fadd(u,v) = [(u[1]+v[1])%p,(u[2]+v[2])%p,(u[3]+v[3])%p];
fneg(u) = [(-u[1])%p,(-u[2])%p,(-u[3])%p];
fscl(c,u) = [(c*u[1])%p,(c*u[2])%p,(c*u[3])%p];
fmul(u,v) = [(u[1]*v[1]+36*(u[2]*v[3]+u[3]*v[2]))%p,(u[1]*v[2]+u[2]*v[1]+36*u[3]*v[3])%p,(u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p];
finv(u) = {my(a,b,c,n,ni); a=u[1]; b=u[2]; c=u[3]; n=lift(Mod(a^3+36*b^3+36^2*c^3-3*36*a*b*c,p)); if(n==0,error("f")); ni=lift(Mod(n,p)^(-1)); [(a^2-36*b*c)*ni%p,(36*c^2-a*b)*ni%p,(b^2-a*c)*ni%p]};
is_cube(r) = (lift(Mod(r,p)^((p-1)/3))==1);
cube_root_fp3(r) = {my(r1,r2,r3,c); r=r%p; if(r==0,return([0,0,0])); r1=r; if(is_cube(r1),c=lift(Mod(r1,p)^5);return([c,0,0])); r2=(r*lift(Mod(36,p)^(-1)))%p; if(is_cube(r2),c=lift(Mod(r2,p)^5);return([0,c,0])); r3=(r*lift(Mod(6,p)^(-1)))%p; if(is_cube(r3),c=lift(Mod(r3,p)^5);return([0,0,c])); error("x")};

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

jac_cp(aa,bb) = {my(h,cp); h=Mod(1,p)*x^6+Mod(aa,p)*x^3+Mod(bb,p); cp=trap(,0,hyperellcharpoly(h)); cp};

\\ Does the degree-4 charpoly cp equal (T^2-t1*T+p)*(T^2-t2*T+p)?
matches_pair(cp, t1, t2) = {
  my(target); if(cp==0, return(0));
  target = (x^2-t1*x+p)*(x^2-t2*x+p);
  cp == target
};

\\ ================================================================
print("================================================================");
print("Full-charpoly search, p=43, E1: b1=7 fixed. b2 in 1..42.");
print("Distribution check: how many DISTINCT char-polys occur overall?");
print("================================================================");
{
  b1=7; E1=ellinit([0,b1],p); n1=ellcard(E1); t1=p+1-n1;
  print("  E1: b=",b1,"  #E1=",n1,"  trace=",t1);

  \\ First pass: collect all valid (r1,r2,k,aa,bb,cp) once (independent of b2).
  all=List();
  for(r1=1,p-1,
    av=cube_root_fp3(r1);
    for(r2=1,p-1,
      bv=cube_root_fp3(r2);
      for(kk=0,2,
        bvk=if(kk==0,bv,if(kk==1,fscl(z3,bv),fscl(z3sq,bv)));
        s0=fadd(av,bvk); q0=fmul(av,bvk);
        res=trap(,[-1,-1],richelot_gen(s0,q0));
        if(res[1]!=-1,
          cp=jac_cp(res[1],res[2]);
          if(cp!=0, listput(all,[r1,r2,kk,res[1],res[2],cp]))
        )
      )
    )
  );
  print("  total valid nonsingular covers: ",#all);
  distinct = Set(Vec(apply(v->v[6], all)));
  print("  distinct char polys among them: ",#distinct);
  print("");

  found_any = 0;
  for(b2=1,p-1,
    if(b2==b1, next());
    E2=ellinit([0,b2],p); n2=ellcard(E2); t2=p+1-n2;
    nhits=0; hitlist=List();
    for(i=1,#all,
      v=all[i];
      if(matches_pair(v[6], t1, t2), nhits++; listput(hitlist,v));
    );
    print("  b2=",b2,"  t2=",t2,"  STRICT hits (full charpoly match): ",nhits);
    if(nhits>0,
      found_any=1;
      for(i=1,min(#hitlist,5),
        v=hitlist[i];
        print("     r1=",v[1]," r2=",v[2]," k=",v[3],
              "  r1/(-b1)=",lift(Mod(v[1],p)/Mod(-b1,p)),
              "  r2/(-b2)=",lift(Mod(v[2],p)/Mod(-b2,p)));
      )
    );
  );
  if(!found_any, print("  NO b2 in 1..42 produced a strict full-charpoly match."));
}
print("================================================================");
print("Done.");
print("================================================================");

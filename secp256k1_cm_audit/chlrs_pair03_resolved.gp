\\ chlrs_pair03_resolved.gp
\\
\\ Follow-up to 2026-07-27's Test 3 (howe_5pairs_v2.gp), which declared secp256k1
\\ pair (0,3) [E0: y^2=x^3+7, E3: y^2=x^3-7] "BLOCKED: degenerate, d=-1, sv=0,
\\ Delta=0". That is true only for the literal branch d=-1. -1 has THREE cube
\\ roots in F_p when p = 1 mod 3 (true for p_secp): d in {-1, -z3, -z3^2} where
\\ z3 is a primitive cube root of unity. d=-1 gives sv=alpha+beta=(1+d)*alpha=0
\\ (degenerate). The other two roots give (1+d) != 0 and should be non-degenerate.
\\ This script tests those branches directly on the real secp256k1 prime.
\\
\\ Run: gp -q chlrs_pair03_resolved.gp

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
  if(nrm == 0, return(0));
  ni = lift(Mod(nrm, pp)^(-1));
  [((a^2 - rr*b*c)*ni) % pp, ((rr*c^2 - a*b)*ni) % pp, ((b^2 - a*c)*ni) % pp]
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
  if(D0inv == 0, return([-1,-1]));

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

jac_order(aa, bb, pp) = {
  my(hh, cp, nj, dsc);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  dsc = lift(poldisc(hh));
  if(dsc == 0, return(-1));
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  nj
};

pp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
b1 = 7; b2 = lift(Mod(-7,pp));
E1 = ellinit([0,b1],pp); E2 = ellinit([0,b2],pp);
t1 = pp+1-ellcard(E1); t2 = pp+1-ellcard(E2);
print("p_secp pair (0,3): b1=",b1," b2=",b2,"  t1=",t1,"  t2=",t2);
print("t1 == -t2 ? ", t1==(-t2)%pp || t1+t2==0);

z3 = lift(polrootsmod(x^2+x+1,pp)[1]);
print("z3 (primitive cube root of unity) = ",z3);
print("z3^3 mod p = ",lift(Mod(z3,pp)^3));

rr = (-b1)%pp;  \\ alpha^3 = rr = -7
print("Confirm -1 is a cube in F_p_secp: (-1)^((p-1)/3) mod p = ", lift(Mod(-1,pp)^((pp-1)/3)));

d_degenerate = lift(Mod(-1,pp));
d_branch1 = lift(Mod(-1,pp)*Mod(z3,pp));
d_branch2 = lift(Mod(-1,pp)*Mod(z3,pp)^2);
print("Three cube roots of -1 mod p: ", [d_degenerate, d_branch1, d_branch2]);
print("Check cubes: ", [lift(Mod(d_degenerate,pp)^3), lift(Mod(d_branch1,pp)^3), lift(Mod(d_branch2,pp)^3)], " (expect all -1 = ",lift(Mod(-1,pp)),")");

results = [];
{
for(bi=1, 2,
  dbr = if(bi==1, d_branch1, d_branch2);
  sv = [0, (1+dbr)%pp, 0];
  qv = [0, 0, dbr%pp];
  print("");
  print("--- branch ",bi," : d=",dbr," ---");
  res = richelot(sv, qv, z3, rr, pp);
  print("  richelot -> (aa,bb) = ",res);
  results = concat(results, [res]);
  if(res[1] != -1,
    hh = Mod(1,pp)*x^6 + Mod(res[1],pp)*x^3 + Mod(res[2],pp);
    dsc = lift(poldisc(hh));
    print("  smooth (disc != 0)? ", dsc != 0);
    print("  (hyperellcharpoly skipped: overflows at 256-bit p; verified equivalent");
    print("   branch construction numerically at p=43 toy scale in chlrs_forward_map_sweep.gp)");
  ,
    print("  DEGENERATE (Delta=0 or non-F_p cover)");
  )
)
}
{
if(results[1][1]!=-1 && results[2][1]!=-1,
  print("");
  print("branch1 aa = ",results[1][1],"  branch2 aa = ",results[2][1]);
  print("branch1 bb = ",results[1][2],"  branch2 bb = ",results[2][2]);
  print("aa1 == -aa2 mod p (same bb)? ", (results[1][1]+results[2][1])%pp==0 && results[1][2]==results[2][2]);
);
}
print("");
print("Done.");

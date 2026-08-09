\\ chlrs_forward_map_sweep.gp
\\
\\ Thread 3 (CHLRS forward map). Prior sessions (2026-07-27, howe_5pairs_v2.gp
\\ Test 1/2) assumed the forward map is:
\\   E1: y^2=x^3+b1, E2: y^2=x^3+b2=d^3*b1  =>  alpha=cbrt(-b1), beta=d*alpha
\\ fed directly into the Z/3Z Richelot dual (sv=alpha+beta, qv=alpha*beta).
\\ Test 1 (p=43, b1=7,b2=13,d=2) was logged as "CORRECT" because it reproduced
\\ literature constants (a=41,b=5) -- but that check never verified #Jac(cover)
\\ against #E1*#E2. This script closes that gap: brute-force over (b1,d) and
\\ over the 3 possible cube-root branches for beta (d*alpha, d*alpha*z3,
\\ d*alpha*z3^2), and check the resulting cover's Jacobian order against all
\\ 4 sign combinations of (t1,t2) (to allow for twist ambiguity).
\\
\\ Run: gp -q chlrs_forward_map_sweep.gp

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

\\ Given (aa,bb), return #Jac(y^2=x^6+aa*x^3+bb) or -1 if not smooth.
jac_order(aa, bb, pp) = {
  my(hh, cp, nj, dsc);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  dsc = lift(poldisc(hh));
  if(dsc == 0, return(-1));
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  nj
};

\\ ================================================================
\\ Brute force over p=43: all (b1,d) pairs, all 3 beta branches.
\\ ================================================================
pp = 43;
z3 = lift(polrootsmod(x^2+x+1,pp)[1]);
z3sq = lift(Mod(z3,pp)^2);
print("p=",pp,"  z3=",z3);

nmatches = 0;
ntested = 0;
{
for(b1=1, pp-1,
  rr = (-b1)%pp;
  \\ alpha = cube root of rr = -b1, as field generator [0,1,0]
  E1 = ellinit([0,b1],pp);
  t1 = pp+1-ellcard(E1);
  for(d=1, pp-1,
    b2 = lift(Mod(d,pp)^3*b1);
    if(b2==b1, next);  \\ same curve
    E2 = ellinit([0,b2],pp);
    t2 = pp+1-ellcard(E2);
    for(branch=0, 2,
      dbr = lift(Mod(d,pp)*Mod(z3,pp)^branch);
      sv = [0, (1+dbr)%pp, 0];
      qv = [0, 0, dbr%pp];
      ntested++;
      res = richelot(sv, qv, z3, rr, pp);
      if(res[1] != -1,
        nj = jac_order(res[1], res[2], pp);
        if(nj > 0,
          targets = [(pp+1-t1)*(pp+1-t2), (pp+1-t1)*(pp+1+t2),
                     (pp+1+t1)*(pp+1-t2), (pp+1+t1)*(pp+1+t2)];
          for(ti=1, 4,
            if(nj == targets[ti],
              nmatches++;
              if(nmatches <= 20,
                print("  MATCH b1=",b1," b2=",b2," d=",d," branch=",branch,
                      "  (a,b)=",res," t1=",t1," t2=",t2," target#",ti," nj=",nj);
              );
            )
          );
        )
      )
    )
  )
);
}
print("Tested ",ntested," (b1,d,branch) triples over F_",pp,": ",nmatches," Jacobian-order matches.");

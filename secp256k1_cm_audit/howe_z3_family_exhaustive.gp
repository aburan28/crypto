\\ howe_z3_family_exhaustive.gp
\\
\\ Thread 3 (CHLRS Igusa forward map / Howe cover construction).
\\
\\ Prior work (Thread 2, 2026-07-27) used a "Z/3Z Richelot" transform
\\ (richelot() below, lifted from howe_5pairs_v2.gp / howe_richelot_v5.gp)
\\ to try to construct the explicit genus-2 Howe cover of a target pair
\\ (E1, E2). It restricted the Richelot input to a 2-parameter family
\\ sv=[0,s,0], qv=[0,0,q] (s,q in F_p) for a fixed cubic-extension
\\ modulus rr, i.e. G1*G2*G3 = (x^3-r1)(x^3-r2) with alpha=r1^(1/3),
\\ beta=r2^(1/3) restricted to be scalar multiples of one fixed F_{p^3}
\\ generator. Thread 2's p=43 "reference" test (sv=[0,3,0], qv=[0,0,2))
\\ was accepted as CORRECT purely because it matched a previously
\\ hard-coded (a,b)=(41,5) -- i.e. a regression check against itself,
\\ never checked against the true target #Jac = #E1 * #E2.
\\
\\ THIS SCRIPT closes that gap: it exhaustively searches the ENTIRE
\\ (s,q) grid (all s,q in F_p), across several representative choices
\\ of rr spanning both non-cube classes of F_p*, and checks the
\\ resulting Richelot-dual curve's true Jacobian order (via
\\ hyperellcharpoly) against the true target #E1(F_p)*#E2(F_p).
\\
\\ RESULT (p=43, E1: y^2=x^3+7, E2: y^2=x^3+13):
\\   target #E1*#E2 = 1767
\\   the p=43 "reference" curve y^2=x^6+41x^3+5 has #Jac = 2169 != 1767.
\\   Exhaustive search: 0 / 14112 F_p-rational (s,q,rr) hits reach 1767,
\\   across rr in {36,7,2,3,5,10,15,20} (spans both non-cube classes).
\\ CONCLUSION: the restricted "Z/3Z scalar-multiple" 2-parameter slice
\\ of the Richelot transform CANNOT reach the Howe cover of (E1,E2) for
\\ this pair -- not a bug in one test, a structural dead end for the
\\ whole family. See RESEARCH_AUTOLAB_LOG.md 2026-08-08 for the writeup
\\ and next-step proposal (general sv,qv in F_{p^3}, not just scalar
\\ multiples of a fixed generator).
\\
\\ Run: gp -q howe_z3_family_exhaustive.gp

default(parisize, 256000000);

\\ ================================================================
\\ F_{p^3} arithmetic (copied from howe_5pairs_v2.gp, unmodified).
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
  if(nrm == 0, return(-1));
  ni = lift(Mod(nrm, pp)^(-1));
  [((a^2 - rr*b*c)*ni) % pp, ((rr*c^2 - a*b)*ni) % pp, ((b^2 - a*c)*ni) % pp]
};

\\ Z/3Z Richelot dual: G1*G2*G3 = (x^3-r1)(x^3-r2) with r1=alpha^3=sv,qv-derived.
\\ Returns [aa,bb] for the dual curve y^2=x^6+aa*x^3+bb if it descends to F_p,
\\ else [-1,-1].
richelot(sv, qv, z3, rr, pp) = {
  my(z3sq, G1c, G1x, G2c, G2x, G3c, G3x, D0, D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q3r,Q6r,lc,lcinv,aa,bb);
  z3sq = lift(Mod(z3, pp)^2);
  G1c = qv;          G1x = f3neg(sv, pp);
  G2c = f3scl(z3sq, qv, pp);  G2x = f3neg(f3scl(z3, sv, pp), pp);
  G3c = f3scl(z3, qv, pp);   G3x = f3neg(f3scl(z3sq, sv, pp), pp);
  D0 = f3add(f3add(
    f3add(f3mul(G2x,G3c,rr,pp), f3neg(f3mul(G3x,G2c,rr,pp),pp), pp),
    f3neg(f3mul(G1x, f3add(G3c, f3neg(G2c,pp), pp), rr, pp), pp), pp),
    f3mul(G1c, f3add(G3x, f3neg(G2x,pp), pp), rr, pp), pp);
  if(D0[1]==0 && D0[2]==0 && D0[3]==0, return([-1,-1]));
  D0inv = f3inv(D0, rr, pp);
  if(D0inv==-1, return([-1,-1]));
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
  Q3r=f3add(f3add(f3mul(P1,H3x2,rr,pp),f3mul(P2,H3x1,rr,pp),pp),f3mul(P3,H3x0,rr,pp),pp);
  Q6r=f3mul(P4,H3x2,rr,pp);
  if(Q0r[2]!=0||Q0r[3]!=0||Q3r[2]!=0||Q3r[3]!=0||Q6r[2]!=0||Q6r[3]!=0,return([-1,-1]));
  lc = Q6r[1]; if(lc==0, return([-1,-1]));
  lcinv = lift(Mod(lc,pp)^(-1));
  aa = (Q3r[1]*lcinv)%pp; bb=(Q0r[1]*lcinv)%pp;
  [aa, bb]
};

\\ ================================================================
\\ Exhaustive check at p=43: E1: y^2=x^3+7, E2: y^2=x^3+13.
\\ Sweeps ALL (s,q) in F_p x F_p, for 8 representative rr values
\\ spanning both non-cube classes of F_43*.
\\ ================================================================
pp=43; z3=6;
E1=ellinit([0,7],pp); E2=ellinit([0,13],pp);
target=ellcard(E1)*ellcard(E2);
print("p=43: E1: y^2=x^3+7, E2: y^2=x^3+13");
print("target #E1*#E2 = ",target);

nhits=0; nmatch=0; matches=[]; rrlist=[36,7,2,3,5,10,15,20];
gettime();
for(ri=1,#rrlist, rr=rrlist[ri]; for(s=0,pp-1, for(q=0,pp-1, {
  res=richelot([0,s,0],[0,0,q],z3,rr,pp);
  if(res[1]!=-1,
    nhits++;
    hh=Mod(1,pp)*x^6+Mod(res[1],pp)*x^3+Mod(res[2],pp);
    nj=iferr(subst(hyperellcharpoly(hh),x,1), E, -1);
    if(nj==target, nmatch++; matches=concat(matches,[[rr,s,q,res[1],res[2]]]));
  );
})));
print("time ms=",gettime());
print("F_p-rational hits=",nhits," / ",pp^2*#rrlist);
print("matches to true target=",nmatch);
print("match list (rr,s,q,aa,bb): ",matches);
if(nmatch==0, print("NEGATIVE RESULT: the restricted (s,q) scalar-multiple slice of the Z/3Z Richelot family never reaches the Howe cover of (E1,E2), for any of the ",#rrlist," representative rr choices tried."));
print("");
print("Done.");

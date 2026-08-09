\\ ============================================================
\\ Mestre's algorithm, forward and backward: a working PARI port
\\ ============================================================
\\
\\ Companion to RESEARCH_MESTRE_HOWE.md §7. Closes two of the three
\\ "open implementation task" items listed there:
\\
\\   1. "Igusa invariants from curve coefficients" (High value,
\\      medium cost) -- implemented below as `igusa_clebsch(f)`,
\\      via Mestre's own Ueberschiebung (transvectant) recipe
\\      instead of hand-transcribing the Cardona-Quer polynomials.
\\
\\   3. "Mestre's Step 2 (conic + sextic)" (Medium value, high
\\      cost) -- implemented below as `reconstruct_curve(...)`.
\\
\\ Item 2 ("Igusa invariants of (E x E^t)/Gamma_alpha", i.e. the
\\ actual gluing formula) is NOT addressed here -- that is a
\\ separate, harder problem (computing invariants of a target
\\ Jacobian we don't have a curve model for yet). This script only
\\ gives us both ends of the Mestre pipeline once those invariants
\\ are in hand.
\\
\\ SOURCE: both formulas below are direct, verified transcriptions
\\ of SageMath's implementation (GPL-2.0+), not independently
\\ re-derived from the classical papers:
\\   - sage/schemes/hyperelliptic_curves/invariants.py
\\     (functions: diffxy, differential_operator, Ueberschiebung,
\\      ubs, clebsch_to_igusa, igusa_clebsch_invariants)
\\   - sage/schemes/hyperelliptic_curves/mestre.py
\\     (functions: Mestre_conic, HyperellipticCurve_from_invariants)
\\ Original references: Mestre, "Construction de courbes de genre 2
\\ à partir de leurs modules", 1991; Lauter-Yang 2001 pp. 956-957
\\ (cross-referenced with van Wamelen 1999 to fix typos, per the
\\ Sage docstrings).
\\
\\ Validated against Sage's own docstring test vectors:
\\   clebsch_to_igusa(2,3,4,5)             == (-240, 17370, 231120, -103098906)
\\   igusa_clebsch_invariants(x^6+1)       == (-240, 1620, -119880, -46656)
\\   igusa_clebsch_invariants(x^6+x^5+x^4+x^2+2)
\\                                          == (-496, 6220, -955932, -1111784)
\\   Mestre_conic([1,2,3,4])               matches Sage's conic coefficients
\\     exactly after clearing denominators (u^2=-2572155000 etc, verified
\\     by hand in the dev session, not re-asserted below to keep this
\\     script finite-field-only)
\\   Mestre_conic([GF(7)(10),1,2,3])       matches Sage's F_7 conic exactly
\\     (0, -2, -1, -2, 2, -3) for (u^2,uv,v^2,uw,vw,w^2)
\\
\\ NEW round-trip test (not in Sage, since Sage doesn't need one --
\\ its two directions are independently trusted; ours needed
\\ cross-checking against each other): for several toy curves over
\\ F_97 and F_1009, compute Igusa-Clebsch invariants, reconstruct a
\\ curve via Mestre, recompute its invariants, and check they agree
\\ up to the projective scaling Mestre's algorithm is only defined
\\ up to (I2^5/I10, I4^5/I10^2, I6^5/I10^3 are scale-invariant).
\\
\\ Run: gp -q mestre_reconstruction.gp

default(parisize, 256000000);
default(timer, 0);

\\ ============================================================
\\ FORWARD: sextic -> Igusa-Clebsch invariants, via transvectants
\\ ============================================================

\\ Partial derivative of a bivariate form f (in variables a,b),
\\ dxa times w.r.t. a and dxb times w.r.t. b.
partial(f, dxa, dxb) = {
  my(h = f);
  for(i=1, dxa, h = deriv(h, a));
  for(i=1, dxb, h = deriv(h, b));
  h
};

\\ Ueberschiebung / transvectant (f,g)_k of Mestre 1991 p.315:
\\ expand (fx*gy - fy*gx)^k via the binomial theorem and replace
\\ each fx^{k-i}fy^i (resp. gx^i gy^{k-i}) monomial by the actual
\\ mixed partial derivative of f (resp. g) of that order.
ub(f, degf, g, degg, k) = {
  my(Cc, acc);
  Cc = (degf-k)! * (degg-k)! / (degf! * degg!);
  acc = 0;
  for(i=0, k, acc += (-1)^i * binomial(k,i) * partial(f,k-i,i) * partial(g,i,k-i));
  Cc*acc
};

\\ Mestre p.317's named transvectant chain -> Clebsch invariants A,B,C,D.
clebsch_invariants_hom(f) = {
  my(i4, Delta, y1, y2, y3, A, B, Cc, D);
  i4    = ub(f,6,f,6,4);
  Delta = ub(i4,4,i4,4,2);
  y1    = ub(f,6,i4,4,4);
  y2    = ub(i4,4,y1,2,2);
  y3    = ub(i4,4,y2,2,2);
  A = ub(f,6,f,6,6);
  B = ub(i4,4,i4,4,4);
  Cc = ub(i4,4,Delta,4,4);
  D = ub(y3,2,y1,2,2);
  [A,B,Cc,D]
};

clebsch_to_igusa(A,B,Cc,D) = {
  [-120*A,
   -720*A^2 + 6750*B,
   8640*A^3 - 108000*A*B + 202500*Cc,
   -62208*A^5 + 972000*A^3*B + 1620000*A^2*Cc - 3037500*A*B^2 - 6075000*B*Cc - 4556250*D]
};

\\ f: bivariate homogeneous degree-6 form in a,b (or pass a dehomogenized
\\ univariate poly in x via homogenize6 first). Returns (I2,I4,I6,I10).
igusa_clebsch(f) = { my(v = clebsch_invariants_hom(f)); clebsch_to_igusa(v[1],v[2],v[3],v[4]) };

\\ homogenize a univariate poly h(x) of degree <= 6 into a bivariate form in a,b
homogenize6(h) = {
  my(v = vector(7, i, polcoeff(h,i-1)));
  sum(i=0, 6, v[i+1]*a^i*b^(6-i))
};

\\ ============================================================
\\ BACKWARD: Igusa-Clebsch invariants -> genus-2 curve (Mestre)
\\ ============================================================

\\ Mestre's conic matrix L(I2,I4,I6,I10): the genus-2 curve exists
\\ over the base field iff the conic v^T L v = 0 has a rational point.
mestre_L(I2,I4,I6,I10) = {
  my(xx,yy,zz,L);
  xx = 8*(1 + 20*I4/I2^2)/225;
  yy = 16*(1 + 80*I4/I2^2 - 600*I6/I2^3)/3375;
  zz = -64*(-10800000*I10/I2^5 - 9 - 700*I4/I2^2 + 3600*I6/I2^3
           + 12400*I4^2/I2^4 - 48000*I4*I6/I2^5)/253125;
  L = [xx+6*yy, 6*xx^2+2*yy, 2*zz;
       6*xx^2+2*yy, 2*zz, 9*xx^3+4*xx*yy+6*yy^2;
       2*zz, 9*xx^3+4*xx*yy+6*yy^2, 6*xx^2*yy+2*yy^2+3*xx*zz];
  [xx,yy,zz,L]
};

bil(L,v,w) = v~ * L * w;
quad(L,v) = bil(L,v,v);

\\ Find a point on the conic v^T L v = 0 over F_p in O(p): dehomogenize
\\ at w=1 and complete the square in v for each u.
find_conic_point(L, p) = {
  my(a,b,c,disc,sq,v);
  for(u=0, p-1,
    a = L[2,2]; b = 2*L[1,2]*u+2*L[2,3]; c = L[1,1]*u^2+2*L[1,3]*u+L[3,3];
    if(a != 0,
      disc = b^2-4*a*c;
      if(issquare(disc, &sq),
        v = (-b+sq)/(2*a);
        return(Mod([u,lift(v),1]~,p))
      )
    ,
      if(b != 0,
        v = -c/b;
        return(Mod([u,lift(v),1]~,p))
      )
    )
  );
  \\ fallback: full brute force including the line at infinity w=0 (rare)
  for(u0=0, p-1, for(u1=0, p-1, for(u2=0, p-1,
    if(u0==0 && u1==0 && u2==0, next);
    v = Mod([u0,u1,u2]~, p);
    if(quad(L,v) == 0, return(v))
  )));
  error("no conic point found over F_", p)
};

\\ Parametrize the conic through P0: for direction D = s*e1 + e2 (e1,e2 a
\\ basis of a plane not containing P0), the second intersection of the
\\ line (P0,D) with the conic is X(s) = (D^T L D)*P0 - 2*(P0^T L D)*D.
\\ As s ranges over F_p (plus s=infinity), X traces the whole conic.
mestre_parametrize(L, P0, p) = {
  my(cands = [Mod([1,0,0]~,p), Mod([0,1,0]~,p), Mod([0,0,1]~,p)], basis = [], s = 's);
  for(i=1,3,
    my(c = cands[i], indep = 0);
    for(j=1,3, for(k=j+1,3, if(P0[j]*c[k] - P0[k]*c[j] != 0, indep=1)));
    if(indep, basis = concat(basis, [c]));
  );
  my(e1 = basis[1], e2 = basis[2], D = s*e1 + e2);
  my(c1 = quad(L, D));
  my(c2 = 2*bil(L, P0, D));
  my(X = c1*P0 - c2*D);
  [lift(X[1]), lift(X[2]), lift(X[3])]
};

\\ Mestre's c_ijk sextic reconstruction (Mestre 1991 p.321/332, via
\\ Lauter-Yang p.957) from x,y,z and the conic parametrization F1,F2,F3.
mestre_sextic(xx,yy,zz,F1,F2,F3) = {
  my(c111,c112,c113,c122,c123,c133,c222,c223,c233,c333);
  c111 = 12*xx*yy - 2*yy/3 - 4*zz;
  c112 = -18*xx^3 - 12*xx*yy - 36*yy^2 - 2*zz;
  c113 = -9*xx^3 - 36*xx^2*yy - 4*xx*yy - 6*xx*zz - 18*yy^2;
  c122 = c113;
  c123 = -54*xx^4 - 36*xx^2*yy - 36*xx*yy^2 - 6*xx*zz - 4*yy^2 - 24*yy*zz;
  c133 = -27*xx^4/2 - 72*xx^3*yy - 6*xx^2*yy - 9*xx^2*zz - 39*xx*yy^2 - 36*yy^3 - 2*yy*zz;
  c222 = -27*xx^4 - 18*xx^2*yy - 6*xx*yy^2 - 8*yy^2/3 + 2*yy*zz;
  c223 = 9*xx^3*yy - 27*xx^2*zz + 6*xx*yy^2 + 18*yy^3 - 8*yy*zz;
  c233 = -81*xx^5/2 - 27*xx^3*yy - 9*xx^2*yy^2 - 4*xx*yy^2 + 3*xx*yy*zz - 6*zz^2;
  c333 = 27*xx^4*yy/2 - 27*xx^3*zz/2 + 9*xx^2*yy^2 + 3*xx*yy^3 - 6*xx*yy*zz + 4*yy^3/3 - 10*yy^2*zz;
  c111*F1^3 + c112*F1^2*F2 + c113*F1^2*F3 + c122*F1*F2^2 + c123*F1*F2*F3
    + c133*F1*F3^2 + c222*F2^3 + c223*F2^2*F3 + c233*F2*F3^2 + c333*F3^3
};

\\ Full backward pass: (I2,I4,I6,I10) over F_p -> univariate sextic h1(x)
\\ over F_p with those Igusa-Clebsch invariants (up to scaling).
reconstruct_curve(I2,I4,I6,I10,p) = {
  my(xyzL = mestre_L(I2,I4,I6,I10));
  my(xx=xyzL[1], yy=xyzL[2], zz=xyzL[3], L=Mod(xyzL[4],p));
  my(P0 = find_conic_point(L, p));
  my(F = mestre_parametrize(L, P0, p));
  my(raw = mestre_sextic(Mod(xx,p),Mod(yy,p),Mod(zz,p), F[1],F[2],F[3]));
  \\ `raw`'s PARI "main variable" can be a *different*, higher-priority
  \\ variable than 's' (e.g. 'a', from the forward-direction functions
  \\ above) even though raw's actual content only depends on s -- PARI
  \\ assigns variable priority by first use in the session, not by which
  \\ variables a polynomial actually depends on, and poldegree/polcoeff
  \\ silently read the wrong axis (degree 0) when that happens. Force
  \\ extraction w.r.t. s explicitly to sidestep this.
  my(clean = sum(k=0, 6, polcoeff(raw,k,s)*s^k));
  lift(clean)
};

\\ ============================================================
\\ Round-trip validation: sextic -> invariants -> Mestre -> sextic
\\ ============================================================

ratios(IC) = [IC[1]^5/IC[4], IC[2]^5/IC[4]^2, IC[3]^5/IC[4]^3];

roundtrip_test(h0, p, label) = {
  print("---- ", label, " (p=", p, ") ----");
  print("h0(x) = ", h0);
  my(d0 = lift(Mod(poldisc(h0),p)));
  print("disc(h0) mod p = ", d0, if(d0==0," [DEGENERATE, skipping]",""));
  if(d0 == 0, return(0));

  my(f0hom = Mod(homogenize6(h0), p));
  my(IC0 = igusa_clebsch(f0hom));
  print("Igusa-Clebsch(h0)  = ", lift(IC0));
  if(lift(IC0[4]) == 0, print("  I10=0, skipping (Mestre needs I10 != 0)"); return(0));

  my(h1s = reconstruct_curve(lift(IC0[1]),lift(IC0[2]),lift(IC0[3]),lift(IC0[4]), p));
  print("Reconstructed h1(s) = ", h1s);
  my(f1hom = Mod(homogenize6(h1s), p));
  my(IC1 = igusa_clebsch(f1hom));
  print("Igusa-Clebsch(h1)  = ", lift(IC1));
  if(lift(IC1[4]) == 0,
    print("  Reconstructed curve is DEGENERATE (I10=0) for this basis choice");
    print("  in mestre_parametrize -- documented failure mode, see log.");
    print("");
    return(0)
  );

  my(r0 = ratios(IC0), r1 = ratios(IC1));
  my(ok = r0 == r1);
  print("Scale-invariant ratios match: ", ok);
  print("");
  ok
};

print("================================================================");
print("Sage docstring test vectors (forward and backward, exact match)");
print("================================================================");
print("");
print("clebsch_to_igusa(2,3,4,5)                  = ", clebsch_to_igusa(2,3,4,5));
print("  expected                                 = [-240, 17370, 231120, -103098906]");
print("igusa_clebsch(x^6+1)                       = ", igusa_clebsch(a^6+b^6));
print("  expected                                 = [-240, 1620, -119880, -46656]");
print("igusa_clebsch(x^6+x^5+x^4+x^2+2)           = ", igusa_clebsch(a^6+a^5*b+a^4*b^2+a^2*b^4+2*b^6));
print("  expected                                 = [-496, 6220, -955932, -1111784]");
print("");

print("================================================================");
print("Mestre round-trip validation (new -- not in Sage's own test suite)");
print("================================================================");
print("");

results = [];
results = concat(results, roundtrip_test(x^6+3*x^5+5*x^4+7*x^3+11*x^2+13*x+17, 97,   "toy A"));
results = concat(results, roundtrip_test(x^6+x^5+x^4+x^3+x^2+x+1,             97,   "toy B"));
results = concat(results, roundtrip_test(x^6-2*x^3+4*x-1,                     1009, "toy C, larger p"));
results = concat(results, roundtrip_test(x^6+3*x^5+5*x^4+7*x^3+11*x^2+13*x+17, 1009, "toy A, larger p"));

print("SUMMARY: ", #results, " trials, ", sum(i=1,#results,results[i]), " matched exactly.");
print("");

\\ ============================================================
\\ Application: forward-direction invariants for the curves that
\\ mestre_scaffold.gp left half-done (it only computed J_10 via
\\ poldisc, noting I_2, I_4, I_6 were "pages-long" and out of scope)
\\ ============================================================
print("================================================================");
print("Full Igusa-Clebsch invariants for mestre_scaffold.gp's curves");
print("================================================================");
print("");

p_toy = 1009;
h_toy = (x^3 + 11)*(x^3 + 515);
IC_toy = igusa_clebsch(Mod(homogenize6(h_toy), p_toy));
print("h_toy = (x^3+11)(x^3+515) over F_", p_toy);
print("  (I2,I4,I6,I10) = ", lift(IC_toy));
print("");

p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
h_secp = (x^3 + 7)*(x^3 + 189);
IC_secp = igusa_clebsch(Mod(homogenize6(h_secp), p_secp));
print("h_secp = (x^3+7)(x^3+189) over F_p_secp: this is a perfectly smooth");
print("genus-2 curve in its own right (I10 below is nonzero, confirming");
print("mestre_scaffold.gp's earlier poldisc check). The degeneracy found in");
print("RESEARCH_AUTOLAB_LOG.md 2026-07-27 is about a DIFFERENT object: the");
print("Richelot correspondence used to try to construct h_secp AS the");
print("Howe-glued cover of E1 x E2 has vanishing discriminant for the pair");
print("(0,3) (d=-1). h_secp existing and being smooth does not by itself");
print("mean it IS that Howe cover -- that identification is still open.");
print("  (I2,I4,I6,I10) = ", lift(IC_secp));
print("  I10 != 0 (h_secp is smooth genus 2, as independently known): ", lift(IC_secp[4]) != 0);

\\ =====================================================================
\\ chlrs_forward_map.gp  --  autolab 2026-08-01, Thread 2 (CHLRS forward map)
\\
\\ Closes the "missing forward map" gap recorded in RESEARCH_AUTOLAB_LOG.md
\\ (2026-07-27 entry, Findings #3 and #7): given two elliptic curves E1, E2
\\ over F_p and an F_p-Galois-equivariant isomorphism psi : E1[2] -> E2[2],
\\ write down the genus-2 curve C with Jac(C) = (E1 x E2)/Gamma_psi
\\ EXPLICITLY as y^2 = h(x).
\\
\\ The CHLRS Igusa-invariant machinery + Mestre reconstruction
\\ (RESEARCH_MESTRE_HOWE.md section 2) is NOT needed: the 2-torsion gluing
\\ admits a direct elementary construction, derived below and validated
\\ against hyperellcharpoly.
\\
\\ ---------------------------------------------------------------------
\\ DERIVATION
\\
\\ A genus-2 curve with two degree-2 maps to elliptic curves (a (2,2)-split
\\ Jacobian) is bielliptic; normalising the extra involution to x -> -x,
\\
\\     C : y^2 = c*(x^2-al1)*(x^2-al2)*(x^2-al3)
\\
\\ and the two elliptic quotients, from (u,v) = (x^2, y) and
\\ (U,W) = (1/x^2, y/x^3), are
\\
\\     E1 : v^2 = c * prod_i (u - al_i)
\\     E2 : W^2 = -c*al1*al2*al3 * prod_i (U - 1/al_i)
\\
\\ So the 2-torsion x-coordinates of the quotients are {al_i} and {1/al_i},
\\ up to affine changes of coordinate on u and U.
\\
\\ Given E1 : y^2 = prod_i (x-a_i), E2 : y^2 = prod_i (x-b_i) and a pairing
\\ a_i <-> b_i, we need affine A(x) = m*x+nu and A'(x) = m'*x+nu' with
\\
\\     al_i = A(a_i)      and     1/al_i = A'(b_i),
\\
\\ i.e. the unique Mobius transformation N with N(a_i) = b_i must factor as
\\ N = A'^{-1} o inv o A. Writing N = [p0,q0; r0,s0] this is possible IFF
\\ r0 != 0, i.e.
\\
\\     *** the gluing is non-degenerate  <=>  N is not affine ***
\\
\\ and then, solving for the parameters,
\\
\\     sh = s0/r0,  Na = prod_i (a_i+sh),  det0 = p0*s0-q0*r0,
\\     t  = Na*det0        (only the square class of t matters),
\\     al_i = t*(a_i+sh),  c = t.
\\
\\ c = m = t is forced mod squares by requiring the quotients to be E1, E2
\\ themselves and not their quadratic twists:
\\     quotient-1 model = c*m^3*prod(x-a_i)              -> needs c*m square
\\     quotient-2 model = -c*al1al2al3*m'^3*prod(x-b_i)  -> needs -c*Na*R square
\\
\\ If r0 = 0 then psi is induced by an isomorphism E1 -> E2 (a scaling of the
\\ 2-torsion x-coordinates), (E1 x E2)/Gamma_psi is DECOMPOSABLE, and no
\\ genus-2 curve exists. This non-degeneracy condition is NOT captured by
\\ (H1)/(H2)/(H3) as implemented in sextic_twist_howe_check.gp.
\\ =====================================================================

default(parisize, 800000000);

\\ --- Mobius matrix sending a[1]->0, a[2]->1, a[3]->oo -----------------
mobius(a) = [ a[2]-a[3], -a[1]*(a[2]-a[3]) ; a[2]-a[1], -a[3]*(a[2]-a[1]) ];

\\ --- the forward map --------------------------------------------------
\\ a, b : paired 2-torsion x-coords (entries t_INTMOD or t_FFELT)
\\ returns 0 if degenerate, else [c, [al1,al2,al3], sh]
glue2(a, b) =
{
  my(M1, M2, N, p0, q0, r0, s0, sh, Na, det0, t, al);
  M1 = mobius(a); M2 = mobius(b);
  N  = M2^(-1) * M1;
  p0 = N[1,1]; q0 = N[1,2]; r0 = N[2,1]; s0 = N[2,2];
  if(r0 == 0, return(0));
  sh   = s0/r0;
  Na   = (a[1]+sh)*(a[2]+sh)*(a[3]+sh);
  det0 = p0*s0 - q0*r0;
  t    = Na*det0;
  al   = [t*(a[1]+sh), t*(a[2]+sh), t*(a[3]+sh)];
  [t, al, sh];
}

sextic(g) = g[1]*(x^2-g[2][1])*(x^2-g[2][2])*(x^2-g[2][3]);

\\ the 6 bijections {1,2,3} -> {1,2,3}
PERMS = [[1,2,3],[2,3,1],[3,1,2],[1,3,2],[3,2,1],[2,1,3]];

\\ trace of y^2 = prod(x-r_i) over F_p
tr_roots(rts, p) =
{
  my(f = (x-rts[1])*(x-rts[2])*(x-rts[3]), c = Vec(lift(f)));
  ellap(ellinit([0, c[2], 0, c[3], c[4]], p));
}

\\ cubic character index of z in F_p*: e with z^((p-1)/3) = z3^e
cubidx(z, p, z3) =
{
  my(v = Mod(z, p)^((p-1)/3));
  if(v == 1, 0, if(v == Mod(z3, p), 1, 2));
}

\\ y^2 = cc*prod(u-al_i) with j = 0: return B with the curve F_p-isomorphic
\\ to y^2 = x^3 + B, or 0 if j != 0.  (al_i may live in F_{p^3}; the
\\ elementary symmetric functions are assumed F_p-rational.)
jzero_class(cc, al, p) =
{
  my(A2, A4, A6, Bp, lin);
  A2 = -(al[1]+al[2]+al[3]);
  A4 =  al[1]*al[2] + al[1]*al[3] + al[2]*al[3];
  A6 = -al[1]*al[2]*al[3];
  A2 = Mod(fpelt(A2, p), p); A4 = Mod(fpelt(A4, p), p); A6 = Mod(fpelt(A6, p), p);
  cc = Mod(fpelt(cc, p), p);
  lin = A4 - A2^2/3;
  if(lin != 0, return([0, lin]));
  Bp = A6 - A2*A4/3 + 2*A2^3/27;
  [lift(cc^3 * Bp), 0];
}

\\ pull an F_p element out of an FFELT (or pass through a t_INTMOD/t_INT)
fpelt(z, p) = if(type(z) == "t_FFELT", polcoef(z.pol, 0), lift(z));

\\ are two j=0 curves y^2=x^3+B, y^2=x^3+B' F_p-isomorphic?  <=> B/B' in (F_p*)^6
same6(B, Bp, p) = Mod(B,p) != 0 && Mod(Bp,p) != 0 && (Mod(B,p)/Mod(Bp,p))^((p-1)/6) == 1;

\\ =====================================================================
print("=====================================================================");
print("SECTION 1 -- forward map validated on generic split curves over F_p");
print("=====================================================================");

{
  my(p, ok, deg, bad, tried, a, b, g, f, cp, t1, t2, tgt);
  p = 1009;
  setrand(1);
  ok = 0; deg = 0; bad = 0; tried = 0;
  for(trial = 1, 300,
    a = vector(3, i, Mod(random(p), p));
    b = vector(3, i, Mod(random(p), p));
    if(#Set(a) < 3 || #Set(b) < 3, next);
    tried++;
    g = glue2(a, b);
    if(g == 0, deg++; next);
    f = sextic(g);
    if(poldisc(f) == 0, deg++; next);
    cp  = hyperellcharpoly(f);
    t1  = tr_roots(a, p);
    t2  = tr_roots(b, p);
    tgt = (x^2 - t1*x + p)*(x^2 - t2*x + p);
    if(cp == tgt, ok++, bad++;
       if(bad <= 3, print("  MISMATCH got ", cp); print("          want ", tgt)))
  );
  print("random split pairs over F_", p, ":  tried=", tried,
        "  charpoly-EXACT=", ok, "  degenerate=", deg, "  mismatch=", bad);
  print("(charpoly-EXACT means hyperellcharpoly(C) = charpoly(E1)*charpoly(E2)");
  print(" on the nose, so the quadratic twist of each quotient is right too.)");
}

\\ =====================================================================
print("");
print("=====================================================================");
print("SECTION 2 -- j=0 sextic twists: which pairs admit an F_p cover?");
print("=====================================================================");
print("");
print("E_k : y^2 = x^3 + B*h^k,  h = primitive 6th root of unity mod p.");
print("e_k in {0,1,2} defined by (-B*h^k)^((p-1)/3) = z3^{e_k}; Frobenius");
print("acts on the roots of x^3+B*h^k by multiplication by z3^{e_k}, and");
print("e_k = 0 <=> x^3+B*h^k splits over F_p <=> E_k[2] is F_p-rational.");
print("");
print("PREDICTION: pair (k1,k2) admits an F_p-rational NON-DEGENERATE");
print("2-torsion gluing  <=>  e_{k1} + e_{k2} = 0 (mod 3).");
print("Reason: all root sets are scalar multiples of each other, so the");
print("even (cyclic) pairings are exactly the isomorphism-induced ones and");
print("give r0 = 0; the odd pairings are non-degenerate but are Galois-");
print("equivariant only when the two Frobenius rotations e_k1, e_k2 are");
print("inverse.  e_k1 = e_k2 != 0 (every quadratic-twist pair) is DEAD.");
print("");

\\ which (k1,k2) does a degree-4 Frobenius charpoly correspond to?
which_pair(cp, tr, p) =
{
  my(res = "none");
  for(i = 1, 6, for(j = i+1, 6,
      if(cp == (x^2 - tr[i]*x + p)*(x^2 - tr[j]*x + p),
         res = Str("(", i-1, ",", j-1, ")"))));
  res;
}

scan_toy(p, B) =
{
  my(z3, h, Bk, tr, ee, w, rootsA, rootsB, g, f, co, cp, tgt, actual,
     nrat, ndeg, nok, nbad, predicted, found);
  z3 = lift(Mod(znprimroot(p), p)^((p-1)/3));
  h  = lift(Mod(znprimroot(p), p)^((p-1)/6));
  Bk = vector(6, k, lift(Mod(B,p) * Mod(h,p)^(k-1)));
  tr = vector(6, k, ellap(ellinit([0,0,0,0,Bk[k]], p)));
  ee = vector(6, k, cubidx(-Bk[k], p, z3));
  print("p = ", p, "   B = ", B, "   h = ", h, "   z3 = ", z3);
  print("  traces t_k = ", tr);
  print("  e_k        = ", ee);
  predicted = []; found = [];
  for(i = 1, 6, for(j = i+1, 6,
      if((ee[i] + ee[j]) % 3 == 0, predicted = concat(predicted, [[i-1,j-1]]))));
  w = ffgen([p,3], 'w);
  for(i = 1, 6, for(j = i+1, 6,
      nrat = 0; ndeg = 0; nok = 0; nbad = 0;
      rootsA = polrootsmod(x^3 + Bk[i]*w^0);
      rootsB = polrootsmod(x^3 + Bk[j]*w^0);
      for(pi = 1, 6,
        g = glue2([rootsA[1], rootsA[2], rootsA[3]],
                  [rootsB[PERMS[pi][1]], rootsB[PERMS[pi][2]], rootsB[PERMS[pi][3]]]);
        if(g == 0, ndeg++; next);
        f  = sextic(g);
        co = Vec(f);
        if(sum(u = 1, #co, if(co[u]^p == co[u], 0, 1)) > 0, next);
        nrat++;
        f = Pol(vector(#co, u, polcoef(co[u].pol, 0))) * Mod(1,p);
        if(poldisc(f) == 0, ndeg++; next);
        cp  = hyperellcharpoly(f);
        tgt = (x^2 - tr[i]*x + p)*(x^2 - tr[j]*x + p);
        if(cp == tgt, nok++, nbad++; actual = which_pair(cp, tr, p))
      );
      if(nok > 0, found = concat(found, [[i-1,j-1]]));
      if(nok > 0 || nrat > 0 || nbad > 0,
         print("  pair (", i-1, ",", j-1, ")  e=(", ee[i], ",", ee[j],
               ")  Fp-rational=", nrat, "  charpoly-EXACT=", nok,
               "  mismatch=", nbad,
               if(nbad > 0, Str("  [curve actually covers pair ", actual, "]"), "")))));
  print("  predicted : ", predicted);
  print("  observed  : ", found);
  print("  MATCH: ", if(predicted == found, "YES", "*** NO ***"));
  print("");
  [predicted, found];
}

{
  my(agree, tot, r);
  agree = 0; tot = 0;
  foreach([[43,7],[67,7],[79,7],[103,7],[127,7],[151,7],[199,7],[211,7],
           [283,7],[307,7],[331,5],[367,3],[439,7],[463,11],[487,7]], pb,
    r = scan_toy(pb[1], pb[2]);
    tot++; if(r[1] == r[2], agree++));
  print("SECTION 2 SUMMARY: criterion matched observation on ", agree,
        " / ", tot, " toy primes");
}

\\ =====================================================================
print("");
print("=====================================================================");
print("SECTION 3 -- secp256k1: the 15 sextic-twist pairs");
print("=====================================================================");

PSEC = 2^256 - 2^32 - 977;
Z3   = lift((-1 + Mod(-3,PSEC)^((PSEC+1)/4)) / 2);
HSEC = lift(-Mod(Z3,PSEC));
BK   = vector(6, k, lift(Mod(7,PSEC) * Mod(HSEC,PSEC)^(k-1)));
EE   = vector(6, k, cubidx(-BK[k], PSEC, Z3));

{
  print("p  = ", PSEC);
  print("z3 = ", Z3);
  print("     ord(z3)=3 ? ", Mod(Z3,PSEC)^3 == 1 && Mod(Z3,PSEC) != 1);
  print("h  = ", HSEC);
  print("     ord(h)=6 ? ", Mod(HSEC,PSEC)^6 == 1 && Mod(HSEC,PSEC)^3 != 1
                            && Mod(HSEC,PSEC)^2 != 1);
  print("");
  print("e_k (cubic character index of -7*h^k):");
  for(k = 1, 6, print("  k=", k-1, "  e_k=", EE[k], "   x^3+7h^k ",
        if(EE[k]==0, "SPLITS over F_p", "irreducible over F_p")));
  print("");
  print("pair   e1 e2  sum   Fp cover?   note");
  for(i = 1, 6, for(j = i+1, 6,
      my(s = (EE[i]+EE[j]) % 3, note = "");
      if(j == i+3, note = "quadratic-twist pair");
      if(EE[i] == 0 && EE[j] == 0, note = "both 2-torsion rational");
      print("(", i-1, ",", j-1, ")    ", EE[i], "  ", EE[j], "   ", s,
            "     ", if(s == 0, "YES       ", "no        "), note)));
}

print("");
print("--- explicit covers for the qualifying secp256k1 pairs ---");

\\ build the explicit sextic for one secp256k1 pair; verify by recovering
\\ the two elliptic quotients and checking they are F_p-isomorphic to
\\ y^2 = x^3+B1 and y^2 = x^3+B2 (a genuine check on the whole F_{p^3}
\\ construction + descent, independent of the derivation).
build_secp(p, B1, B2, e1, e2) =
{
  my(w, rootsA, rootsB, g, f, co, al, cc, q1, q2, res);
  res = [];
  if(e1 == 0 && e2 == 0,
     rootsA = polrootsmod(x^3 + B1, p);
     rootsB = polrootsmod(x^3 + B2, p)
  ,
     w = ffgen(Mod(1,p)*(x^3 + B1), 'w);
     rootsA = polrootsmod(x^3 + Mod(B1,p)*w^0);
     rootsB = polrootsmod(x^3 + Mod(B2,p)*w^0));
  for(pi = 1, 6,
    g = glue2([rootsA[1], rootsA[2], rootsA[3]],
              [rootsB[PERMS[pi][1]], rootsB[PERMS[pi][2]], rootsB[PERMS[pi][3]]]);
    if(g == 0, next);
    f  = sextic(g);
    co = Vec(f);
    if(sum(u = 1, #co, if(fpelt(co[u]^p, p) != fpelt(co[u], p), 1, 0)) > 0, next);
    f  = Pol(vector(#co, u, fpelt(co[u], p)));
    if(Mod(poldisc(f), p) == 0, next);
    al = g[2]; cc = g[1];
    q1 = jzero_class(cc, al, p);
    q2 = jzero_class(-cc*al[1]*al[2]*al[3], [1/al[1], 1/al[2], 1/al[3]], p);
    res = concat(res, [[pi, f, q1, q2]]);
  );
  res;
}

{
  my(r, f, q1, q2, nqual);
  nqual = 0;
  for(i = 1, 6, for(j = i+1, 6,
      if((EE[i]+EE[j]) % 3 != 0, next);
      nqual++;
      print("");
      print("### pair (", i-1, ",", j-1, ")   B1 = 7*h^", i-1, "  B2 = 7*h^", j-1);
      r = build_secp(PSEC, BK[i], BK[j], EE[i], EE[j]);
      if(#r == 0, print("  no F_p-rational non-degenerate gluing found"); next);
      print("  ", #r, " F_p-rational non-degenerate gluing(s) found");
      f = r[1][2]; q1 = r[1][3]; q2 = r[1][4];
      print("  representative (perm #", r[1][1], "):  C : y^2 = h(x)");
      print("    h(x) = ", f);
      print("  quotient 1 is y^2=x^3+B with B = ", q1[1]);
      print("    j=0 (linear coeff vanishes) : ", q1[2] == 0);
      print("    F_p-isomorphic to E_", i-1, " : ", same6(q1[1], BK[i], PSEC));
      print("  quotient 2 is y^2=x^3+B with B = ", q2[1]);
      print("    j=0 (linear coeff vanishes) : ", q2[2] == 0);
      print("    F_p-isomorphic to E_", j-1, " : ", same6(q2[1], BK[j], PSEC));
  ));
  print("");
  print("qualifying secp256k1 pairs: ", nqual, " / 15");
}

\\ =====================================================================
print("");
print("=====================================================================");
print("SECTION 4 -- CLOSED FORM of the j=0 forward map");
print("=====================================================================");
print("");
print("Specialising the construction to E1: y^2=x^3+B1, E2: y^2=x^3+B2:");
print("the roots of x^3+B1 are a*z3^i, those of x^3+B2 are d*a*z3^{-i} for");
print("the odd pairing, so N(z) = d*a^2/z, i.e. s0 = 0, hence sh = 0 and");
print("");
print("      tau := B1 * (B1*B2)^(1/3)          [in F_p iff B1*B2 is a cube]");
print("      C    : y^2 = tau*x^6 + tau^4*B1");
print("           ~ y^2 = tau*(x^6 + B2/B1)     [x-scaling by a 6th power]");
print("");
print("So the explicit Howe-glued cover of two j=0 curves is the CLASSICAL");
print("bielliptic sextic y^2 = c*(x^6 + r); no Igusa inversion or Mestre");
print("reconstruction is involved at all.");
print("");

\\ closed form: returns the sextic, or 0 if B1*B2 is not a cube in F_p
closed_form(p, B1, B2) =
{
  my(cr, tau);
  if(Mod(B1*B2, p)^((p-1)/3) != 1, return(0));
  cr  = lift(sqrtn(Mod(B1*B2, p), 3));
  tau = lift(Mod(B1,p) * cr);
  Mod(tau,p)*x^6 + Mod(tau,p)^4*Mod(B1,p);
}

{
  my(z3, h, Bk, tr, f, cp, tgt, nok, nbad, nskip, tot);
  nok = 0; nbad = 0; nskip = 0; tot = 0;
  foreach([43,67,79,103,127,151,199,211,283,307,331,367,439,463,487,
           523,547,571,607,619,631,643,691,727,739,751,787,811,823,859], p,
    z3 = lift(Mod(znprimroot(p), p)^((p-1)/3));
    h  = lift(Mod(znprimroot(p), p)^((p-1)/6));
    Bk = vector(6, k, lift(Mod(7,p) * Mod(h,p)^(k-1)));
    tr = vector(6, k, ellap(ellinit([0,0,0,0,Bk[k]], p)));
    for(i = 1, 6, for(j = i+1, 6,
        f = closed_form(p, Bk[i], Bk[j]);
        tot++;
        if(f == 0, nskip++; next);
        if(poldisc(f) == 0, nbad++; next);
        cp  = hyperellcharpoly(f);
        tgt = (x^2 - tr[i]*x + p)*(x^2 - tr[j]*x + p);
        if(cp == tgt, nok++, nbad++;
           if(nbad <= 3, print("  CLOSED-FORM MISMATCH p=", p, " pair (", i-1,
                               ",", j-1, ") got ", factor(cp)))))));
  print("closed form over 30 primes: ", tot, " pairs, ",
        nskip, " skipped (B1*B2 not a cube), ",
        nok, " charpoly-EXACT, ", nbad, " WRONG");
}

print("");
print("--- closed form for secp256k1 ---");
{
  my(f, i, j);
  for(i = 1, 6, for(j = i+1, 6,
      f = closed_form(PSEC, BK[i], BK[j]);
      if(f == 0, next);
      print("pair (", i-1, ",", j-1, "):  y^2 = ", f)));
}

\\ =====================================================================
print("");
print("=====================================================================");
print("SECTION 5 -- the two constructions agree on secp256k1");
print("=====================================================================");
print("");
print("y^2 = A*x^6+C and y^2 = A'*x^6+C' are F_p-isomorphic iff C/C' is a");
print("square and (A/C)/(A'/C') is a 6th power (x -> lam*x, y -> mu*y).");
print("");

iso6(A, C, Ap, Cp, p) =
{
  (Mod(C,p)/Mod(Cp,p))^((p-1)/2) == 1 &&
  ((Mod(A,p)/Mod(C,p))/(Mod(Ap,p)/Mod(Cp,p)))^((p-1)/6) == 1;
}

{
  my(r, f1, f2, A, C, Ap, Cp, nagree, ntot);
  nagree = 0; ntot = 0;
  for(i = 1, 6, for(j = i+1, 6,
      if((EE[i]+EE[j]) % 3 != 0, next);
      r  = build_secp(PSEC, BK[i], BK[j], EE[i], EE[j]);
      f1 = r[1][2];
      f2 = closed_form(PSEC, BK[i], BK[j]);
      A  = polcoef(f1, 6); C  = polcoef(f1, 0);
      Ap = lift(polcoef(f2, 6)); Cp = lift(polcoef(f2, 0));
      ntot++;
      if(iso6(A, C, Ap, Cp, PSEC), nagree++;
         print("pair (", i-1, ",", j-1, "): Mobius construction == closed form  (F_p-isomorphic)"),
         print("pair (", i-1, ",", j-1, "): *** DISAGREE ***"))));
  print("");
  print("agreement: ", nagree, " / ", ntot, " qualifying secp256k1 pairs");
}

print("");
print("done.");
quit

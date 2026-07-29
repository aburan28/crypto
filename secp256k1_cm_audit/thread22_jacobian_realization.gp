\\ thread22_jacobian_realization.gp
\\
\\ Thread 22: Which of the 15 sextic-twist pairs (E_i, E_j) have their product
\\ isogeny class REALIZED by a genus-2 Jacobian over F_p?
\\
\\ MOTIVATION (from 2026-07-21 / Thread 18 and 2026-07-26 / Thread 2):
\\   Thread 18 checked Howe's gluing conditions (H1)(H2)(H3) for all 15 pairs of
\\   secp256k1's j=0 sextic twists and found 5/15 qualify. Thread 2 then showed the
\\   F_p Rosenhain construction is permanently blocked for secp256k1 (cubic-residue
\\   obstruction), so we cannot exhibit the cover curve directly at cryptographic size.
\\   The 2026-07-26 run proposed "Thread 22: Richelot search over small proxy primes".
\\
\\ WHAT THIS SCRIPT DOES (and what it does NOT):
\\   Over small proxy primes p = 1 mod 6 (so all six j=0 sextic twists exist over F_p)
\\   we EXHAUSTIVELY enumerate every genus-2 curve over F_p, compute its Frobenius
\\   char poly, and record which of the 15 products P_i(T)*P_j(T) occur.
\\
\\   This answers the ISOGENY-CLASS question ("does the isogeny class of E_i x E_j
\\   contain a Jacobian?"), which is WEAKER than Howe's gluing question ("is E_i x E_j
\\   itself (2,2)-isogenous to a Jacobian?"). Logically:
\\       H1 & H2 & H3  ==>  class realized.
\\   So a pair that PASSES H1-H3 but whose class is NOT realized would be a
\\   contradiction (i.e. a bug). A pair that FAILS H1-H3 but IS realized merely shows
\\   the Howe conditions are sufficient, not necessary. Both outcomes are informative.
\\
\\ COMPLETENESS OF THE ENUMERATION (this is the part that makes the negatives real):
\\   Every genus-2 curve over F_p, p >= 7, has a model y^2 = f(x) with deg f = 6 and
\\   leading coefficient f6 != 0:
\\     - a deg-6 model with f6 != 0 is already of this form;
\\     - a deg-5 model y^2 = g(x): translate x -> x+c so g(0) != 0 (possible, g has
\\       <= 5 roots and p >= 7), then x -> 1/x, y -> y/x^3 gives y^2 = x^6*g(1/x),
\\       of degree 6 with leading coefficient g(0) != 0.
\\   Given f6 != 0 and p > 3, the translation x -> x - f5/(6*f6) sets f5 = 0.
\\   The remaining transformations x -> u*x, y -> e*y send f6 -> f6*u^6/e^2 and
\\   preserve f5 = 0; the set {u^6 * e^-2} equals the squares of F_p^* (6th powers are
\\   squares, and e^-2 already sweeps all squares), so f6 normalises to {1, nu} with
\\   nu any fixed non-residue.
\\   Hence  { f6*x^6 + f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0 : f6 in {1,nu} }
\\   is a COMPLETE (redundant) set of models: 2*p^5 polynomials, of which the
\\   squarefree ones are exactly the genus-2 curves.
\\
\\ CONVENTION CHECK (verified numerically below):
\\   PARI hyperellcharpoly returns chi(T) = T^4 + a1*T^3 + a2*T^2 + p*a1*T + p^2 with
\\   #C(F_p) = p + 1 + a1. The pair (a1, a2) determines chi, so we index on (a1, a2).
\\   For E_i, E_j of traces t_i, t_j:
\\     P_i*P_j = (T^2 - t_i*T + p)(T^2 - t_j*T + p)
\\             => a1 = -(t_i + t_j),  a2 = 2*p + t_i*t_j.
\\
\\ Run: gp -q thread22_jacobian_realization.gp
\\ Optional: gp -q -D... ; primes are set in PRIMES below.

default(parisize, 2000000000);
default(timer, 0);

PRIMES = [7, 13, 19];

\\ ---------------------------------------------------------------------------
\\ Convention self-check: compare hyperellcharpoly against a direct point count.
\\ ---------------------------------------------------------------------------
convention_check() = {
  my(p, f, ch, a1, cnt, v, ok);
  p = 13;
  f = Mod(1,p)*(x^6 + 2*x^4 + 3*x^3 + x + 5);
  ch = hyperellcharpoly(f);
  a1 = polcoeff(ch, 3);
  cnt = 0;
  for(a = 0, p-1,
    v = lift(subst(lift(f), x, a)) % p;
    if(v == 0, cnt += 1, if(kronecker(v, p) == 1, cnt += 2));
  );
  cnt += 2;  \\ two points at infinity: leading coeff 1 is a square
  ok = (cnt == p + 1 + a1);
  print("convention check: chi = ", ch);
  print("  #C(F_p) counted = ", cnt, "   p+1+a1 = ", p+1+a1, "   -> ",
        if(ok, "OK", "MISMATCH"));
  if(!ok, error("hyperellcharpoly convention mismatch"));
};

\\ ---------------------------------------------------------------------------
\\ Sextic twist data for the j=0 family over F_p.
\\ ---------------------------------------------------------------------------
\\ Returns [bs, ts, ords, pats] with
\\   bs[k]   = b0 * u^(k-1),  u a primitive 6th root of unity  (k = 1..6)
\\   ts[k]   = trace of E_k : y^2 = x^3 + bs[k]
\\   ords[k] = #E_k(F_p)
\\   pats[k] = sorted degree pattern of the factorisation of x^3 + bs[k] over F_p
twist_data(p) = {
  my(g, u, bs, ts, ords, pats, E, fac, pat);
  g = znprimroot(p);
  u = lift(g^((p-1)/6));
  \\ Any nonzero b0 gives the same six twist classes (u^k, k=0..5, is a full set of
  \\ coset representatives for F_p^* / (F_p^*)^6). We use b0 = 1.
  bs = vector(6, k, lift(Mod(u, p)^(k-1)));
  ts = vector(6); ords = vector(6); pats = vector(6);
  for(k = 1, 6,
    E = ellinit([0, 0, 0, 0, bs[k]], p);
    ords[k] = ellcard(E);
    ts[k] = p + 1 - ords[k];
    fac = factormod(Mod(1,p)*(x^3 + bs[k]));
    pat = vecsort(vector(matsize(fac)[1], i,
                         poldegree(fac[i,1]) * fac[i,2]));
    pats[k] = pat;
  );
  [bs, ts, ords, pats];
};

\\ ---------------------------------------------------------------------------
\\ Exhaustive genus-2 enumeration over F_p; returns everything we measured.
\\ ---------------------------------------------------------------------------
\\ AMAX  = floor(4*sqrt(p)) bounds |a1|;  BMAX = 6*p bounds |a2|.
\\ hitcnt[a1+AMAX+1][a2+BMAX+1] = number of curve models with that (a1,a2).
enumerate(p) = {
  my(AMAX, BMAX, hitcnt, nu, ncurves, nmodels, g, ch, a1, a2);
  AMAX = floor(4 * sqrt(p));
  BMAX = 6 * p;
  hitcnt = matrix(2*AMAX+1, 2*BMAX+1);
  nu = 0;
  for(c = 2, p-1, if(kronecker(c, p) == -1, nu = c; break));
  ncurves = 0; nmodels = 0;
  for(i6 = 1, 2,
    my(c6 = if(i6 == 1, 1, nu));
    for(c4 = 0, p-1, for(c3 = 0, p-1, for(c2 = 0, p-1, for(c1 = 0, p-1, for(c0 = 0, p-1,
      nmodels++;
      g = Mod(1,p) * (c6*x^6 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0);
      if(poldisc(g) != 0,
        ncurves++;
        ch = hyperellcharpoly(g);
        a1 = polcoeff(ch, 3);
        a2 = polcoeff(ch, 2);
        hitcnt[a1+AMAX+1, a2+BMAX+1]++;
      );
    ))))));
  [hitcnt, AMAX, BMAX, ncurves, nmodels, nu];
};

\\ ---------------------------------------------------------------------------
\\ Context table: the SAME enumeration answers the question for EVERY split
\\ ordinary class E_{t1} x E_{t2}, not just the 15 twist pairs. List the split
\\ classes that no Jacobian realizes -- the exceptions to "split => Jacobian".
\\ ---------------------------------------------------------------------------
\\ Which traces occur over F_p? By Deuring/Waterhouse, for p prime every t with
\\ |t| <= 2*sqrt(p) and p !| t occurs; t = 0 occurs iff p = 3 mod 4. Since
\\ 2*sqrt(p) < p for p >= 5, the only t divisible by p in range is t = 0.
trace_exists(p, t) = if(t != 0, 1, p % 4 == 3);

missing_split_classes(p, hitcnt, AMAX, BMAX) = {
  my(B, nmiss, ntot, aa1, aa2, c);
  B = sqrtint(4*p);
  nmiss = 0; ntot = 0;
  for(t1 = -B, B,
    if(trace_exists(p, t1),
      for(t2 = t1, B,
        if(trace_exists(p, t2),
          aa1 = -(t1 + t2);
          aa2 = 2*p + t1*t2;
          ntot++;
          c = 0;
          if(abs(aa1) <= AMAX && abs(aa2) <= BMAX, c = hitcnt[aa1+AMAX+1, aa2+BMAX+1]);
          if(c == 0,
            nmiss++;
            print("    t1=", t1, ", t2=", t2, "  (a1,a2)=(", aa1, ",", aa2, ")",
                  if(t1 == t2, "   [E x E : square class]", ""));
          );
        );
      );
    );
  );
  print("  -> ", nmiss, " of ", ntot, " split classes over F_", p, " have NO Jacobian.");
};

\\ ---------------------------------------------------------------------------
\\ Main per-prime report.
\\ ---------------------------------------------------------------------------
run_prime(p) = {
  my(td, bs, ts, ords, pats, en, hitcnt, AMAX, BMAX, ncurves, nmodels, nu);
  my(a1, a2, h1, h2, h3, howe, realized, nreal, nhowe, agree, ndistinct);
  my(rows, tag, verdict);

  print("");
  print("================================================================");
  print("p = ", p, "    (p mod 6 = ", p % 6, ")");
  print("================================================================");

  td = twist_data(p);
  bs = td[1]; ts = td[2]; ords = td[3]; pats = td[4];

  print("");
  print("Six j=0 sextic twists  E_k : y^2 = x^3 + b_k   over F_", p);
  print("  k | b_k | #E_k | t_k | 2-torsion pattern of x^3+b_k");
  for(k = 1, 6,
    print("  ", k-1, " | ", bs[k], " | ", ords[k], " | ", ts[k], " | ", pats[k]);
  );
  print("  trace sums (should be 0 for k in {0,2,4} and {1,3,5}): ",
        ts[1]+ts[3]+ts[5], ", ", ts[2]+ts[4]+ts[6]);

  en = enumerate(p);
  hitcnt = en[1]; AMAX = en[2]; BMAX = en[3];
  ncurves = en[4]; nmodels = en[5]; nu = en[6];

  ndistinct = 0;
  for(i = 1, 2*AMAX+1, for(j = 1, 2*BMAX+1,
    if(hitcnt[i,j] > 0, ndistinct++)));

  print("");
  print("Exhaustive genus-2 enumeration over F_", p, " (non-residue nu = ", nu, "):");
  print("  models scanned (f6 in {1,nu}, f5=0, f4..f0 free) = ", nmodels, "  (= 2*p^5)");
  print("  squarefree (= genus-2 curves)                    = ", ncurves);
  print("  distinct Frobenius char polys realized           = ", ndistinct);

  print("");
  print("15 sextic-twist pairs:  Howe(H1,H2,H3) prediction  vs  actual realization");
  print("  pair | H1 | H2 | H3 | Howe | (a1,a2) target      | #models | realized?");
  nreal = 0; nhowe = 0; agree = 1;
  for(i = 1, 6, for(j = i+1, 6,
    h1 = (ords[i] != ords[j]);
    h2 = (pats[i] == pats[j]);
    h3 = (gcd(ords[i], ords[j]) == 1);
    howe = (h1 && h2 && h3);
    a1 = -(ts[i] + ts[j]);
    a2 = 2*p + ts[i]*ts[j];
    my(cnt = 0);
    if(abs(a1) <= AMAX && abs(a2) <= BMAX,
       cnt = hitcnt[a1+AMAX+1, a2+BMAX+1]);
    realized = (cnt > 0);
    if(howe, nhowe++);
    if(realized, nreal++);
    \\ The only logically forbidden cell: Howe says glueable, nothing realizes it.
    if(howe && !realized, agree = 0);
    tag = concat(concat(concat("(", Str(i-1)), concat(",", Str(j-1))), ")");
    print("  ", tag, "  |  ", if(h1,"Y","n"), " |  ", if(h2,"Y","n"),
          " |  ", if(h3,"Y","n"), " |  ", if(howe,"YES"," no "),
          " | (", a1, ", ", a2, ")",
          " | ", cnt, " | ", if(realized, "YES", "no"));
  ));

  \\ -------------------------------------------------------------------------
  \\ Context: the SAME enumeration answers the question for EVERY split ordinary
  \\ class E_{t1} x E_{t2}, not just the 15 twist pairs. Report the full table of
  \\ missing split classes -- these are the exceptions to "split class => Jacobian".
  \\ -------------------------------------------------------------------------
  print("");
  print("All split classes E_{t1} x E_{t2} over F_", p, " that are NOT realized by any Jacobian:");
  missing_split_classes(p, hitcnt, AMAX, BMAX);

  print("");
  print("  pairs passing H1&H2&H3            : ", nhowe, "/15");
  print("  pairs whose class IS realized     : ", nreal, "/15");
  verdict = if(agree,
    "CONSISTENT: every Howe-qualifying pair is realized as a Jacobian isogeny class.",
    "CONTRADICTION: some Howe-qualifying pair has NO Jacobian in its isogeny class.");
  print("  ", verdict);
  [p, nhowe, nreal, agree];
};

\\ ---------------------------------------------------------------------------
print("================================================================");
print("Thread 22: Jacobian realization of j=0 sextic-twist products");
print("Exhaustive genus-2 search over small proxy primes p = 1 mod 6");
print("================================================================");
print("");
convention_check();

SUMMARY = vector(#PRIMES);
for(i = 1, #PRIMES, SUMMARY[i] = run_prime(PRIMES[i]));

print("");
print("================================================================");
print("SUMMARY");
print("================================================================");
print("   p | Howe-qualifying | class-realized | Howe => realized ?");
{
for(i = 1, #SUMMARY,
  my(s = SUMMARY[i]);
  print("  ", s[1], " |      ", s[2], "/15       |     ", s[3], "/15      | ",
        if(s[4], "OK", "VIOLATED"));
);
}
print("");
print("Reading: 'Howe => realized' must always hold (H1&H2&H3 constructs an explicit");
print("(2,2)-glued Jacobian). 'realized' strictly exceeding 'Howe-qualifying' shows the");
print("Howe conditions are SUFFICIENT but NOT NECESSARY for a split isogeny class to");
print("contain a Jacobian -- i.e. counting Howe-glueable pairs UNDERCOUNTS the genus-2");
print("covers available to a cover attack.");
quit();

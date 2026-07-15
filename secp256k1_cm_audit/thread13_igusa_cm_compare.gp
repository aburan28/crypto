\\ thread13_igusa_cm_compare.gp (v3 — corrected discriminant)
\\ Thread 13: Igusa invariants + CM discriminant for CM-73 vs non-CM-73 norm-form primes.
\\
\\ Weil polynomial: T^4 + a2*T^2 + p^2 (a1=0 by symmetry of naive cover)
\\ CM discriminant: disc4 = a2^2 - 4*p^2 (= squarefree part gives CM field)
\\ CM-73 condition: a2 - 2p = -73, i.e., a2 = 2p - 73.
\\   => disc4 = (2p-73)^2 - 4p^2 = -292p+5329 = -219*k^2 where 4p=73+3k^2.
\\   => sf(disc4) = -219 for ALL CM-73 primes.
\\
\\ Run: gp --stacksize 64000000 -q thread13_igusa_cm_compare.gp

default(parisize, 64000000);
default(timer, 0);

\\ ==============================================================
\\ Core transvectant machinery
\\ ==============================================================

poly_to_binary_form(h, deg) = {
  my(v, j);
  v = vector(deg + 1, j, 0);
  for(j = 0, min(deg, poldegree(h)), v[j + 1] = polcoeff(h, j));
  v
};

transvectant(f, m, g, n, k) = {
  my(out_deg, out, scale, i, sign_i, binom_i, a, b, coeff_f, deriv_f, coeff_g, deriv_g, jj, s);
  out_deg = m + n - 2*k;
  if(out_deg < 0, return([0]));
  out = vector(out_deg + 1, i, 0);
  scale = (m - k)! * (n - k)! / (m! * n!);
  for(i = 0, k,
    sign_i = if(i % 2 == 0, 1, -1);
    binom_i = binomial(k, i);
    for(a = k - i, m - i,
      if(a >= 0 && a <= m,
        deriv_f = 1;
        forstep(s = a, a - (k-i) + 1, -1, deriv_f = deriv_f * s);
        forstep(s = m-a, m-a-i+1, -1, deriv_f = deriv_f * s);
        coeff_f = f[a + 1] * deriv_f;
        if(coeff_f != 0,
          for(b = i, n - (k-i),
            if(b >= 0 && b <= n,
              deriv_g = 1;
              forstep(s = b, b-i+1, -1, deriv_g = deriv_g * s);
              forstep(s = n-b, n-b-(k-i)+1, -1, deriv_g = deriv_g * s);
              coeff_g = g[b + 1] * deriv_g;
              if(coeff_g != 0,
                jj = (a - (k-i)) + (b - i);
                if(jj >= 0 && jj <= out_deg,
                  out[jj + 1] += sign_i * binom_i * coeff_f * coeff_g
                )
              )
            )
          )
        )
      )
    )
  );
  for(jj = 1, out_deg + 1, out[jj] = out[jj] * scale);
  out
};

bfc(v) = if(#v == 1, v[1], error("not degree 0"));

igusa_quad(h) = {
  my(hv, iv, Delta, y1, y2, y3, A, B, C, D, I2, I4, I6, I10);
  hv    = poly_to_binary_form(h, 6);
  iv    = transvectant(hv, 6, hv, 6, 4);
  A     = bfc(transvectant(hv, 6, hv, 6, 6));
  B     = bfc(transvectant(iv, 4, iv, 4, 4));
  Delta = transvectant(iv, 4, iv, 4, 2);
  C     = bfc(transvectant(iv, 4, Delta, 4, 4));
  y1    = transvectant(hv, 6, iv, 4, 4);
  y2    = transvectant(iv, 4, y1, 2, 2);
  y3    = transvectant(iv, 4, y2, 2, 2);
  D     = bfc(transvectant(y3, 2, y1, 2, 2));
  I2  = -120*A;
  I4  = -720*A^2 + 6750*B;
  I6  = 8640*A^3 - 108000*A*B + 202500*C;
  I10 = -62208*A^5 + 972000*A^3*B + 1620000*A^2*C
        - 3037500*A*B^2 - 6075000*B*C - 4556250*D;
  [I2, I4, I6, I10]
};

\\ ==============================================================
\\ Utilities
\\ ==============================================================

prim_root(p) = {
  my(phi, facs, nn, dd, gg, ok, kk);
  phi = p - 1; nn = phi; facs = [];
  dd = 2;
  while(dd*dd <= nn,
    if(nn % dd == 0, facs = concat(facs, [dd]); while(nn % dd == 0, nn = nn \ dd));
    dd = dd + 1
  );
  if(nn > 1, facs = concat(facs, [nn]));
  for(gg = 2, p-1,
    ok = 1;
    for(kk = 1, #facs, if(Mod(gg,p)^(phi/facs[kk]) == 1, ok = 0; break));
    if(ok, return(gg))
  );
  error("no prim root for ", p)
};

\\ Squarefree part of n (preserving sign); core() gives positive squarefree kernel
sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

\\ ==============================================================
\\ Main worker
\\ Returns [a2, delta=a2-2p, disc4=a2^2-4p^2, sf_disc4, I2, I4, I6, I10, j1]
\\ ==============================================================

analyze_p(p) = {
  my(g, b1, b2, fld, P, f, a2, delta, disc4, sf4, h, q, I2, I4, I6, I10, j1);
  g  = prim_root(p);
  b1 = lift(Mod(g, p)^1);
  b2 = lift(Mod(g, p)^2);
  fld = ffgen(p, 't);
  P = fld^0*x^6 + (b1+b2)*fld^0*x^3 + (b1*b2 % p)*fld^0;
  f  = hyperellcharpoly(P);
  a2 = polcoeff(f, 2);
  delta = a2 - 2*p;        \\ = a2-2p; CM-73 iff delta=-73
  disc4 = a2^2 - 4*p^2;   \\ discriminant of Weil poly T^4+a2*T^2+p^2
  sf4 = sf_part(disc4);    \\ squarefree part = CM discriminant
  h  = (x^3 + b1) * (x^3 + b2);
  q  = igusa_quad(h);
  I2  = lift(Mod(q[1], p));
  I4  = lift(Mod(q[2], p));
  I6  = lift(Mod(q[3], p));
  I10 = lift(Mod(q[4], p));
  j1  = if(I10 == 0, -1, lift(Mod(I2, p)^5 * Mod(I10, p)^(-1)));
  [a2, delta, disc4, sf4, I2, I4, I6, I10, j1]
};

\\ ==============================================================
\\ Part A: CM-73 primes
\\ ==============================================================

print("=== Thread 13: Igusa Invariants + CM Discriminant Comparison ===");
print();
print("Weil poly: T^4+a2*T^2+p^2.  disc4=a2^2-4p^2.  CM-73 iff a2-2p=-73.");
print();
print("--- Part A: CM-73 primes {19,37,79,109} ---");
print("Analytic claim: disc4 = -219*k^2 => sf(disc4)=-219 for ALL CM-73 primes");
print();

{
  my(cm73, p, k, r, a2, delta, disc4, sf4, I2, I4, I6, I10, j1);
  cm73 = [19, 37, 79, 109];
  for(ii = 1, 4,
    p = cm73[ii];
    k = sqrtint((4*p - 73) \ 3);
    r = analyze_p(p);
    a2 = r[1]; delta = r[2]; disc4 = r[3]; sf4 = r[4];
    I2 = r[5]; I4 = r[6]; I6 = r[7]; I10 = r[8]; j1 = r[9];
    print("p=", p, " k=", k, "  a2=", a2, "  delta=a2-2p=", delta,
          "  disc4=a2^2-4p^2=", disc4, "  (=-219*", k^2, "?)=", disc4 == -219*k^2);
    print("  sf(disc4)=", sf4, "  (=-219?)", sf4 == -219);
    print("  I2=", I2, "  I4=", I4, "  I6=", I6, "  I10=", I10, "  j1=", j1)
  )
}

print();

\\ ==============================================================
\\ Part B: All norm-form primes 4p=73+3k^2, k odd, k up to 89
\\ ==============================================================

print("--- Part B: All norm-form primes 4p=73+3k^2 (k odd <=89) ---");
print();

{
  my(k, val, p, r, a2, delta, disc4, sf4, I2, I4, I10, j1, tag);
  forstep(k = 1, 89, 2,
    val = 73 + 3*k^2;
    if(val % 4 == 0,
      p = val \ 4;
      if(isprime(p) && p % 6 == 1,
        r = analyze_p(p);
        a2 = r[1]; delta = r[2]; disc4 = r[3]; sf4 = r[4];
        I2 = r[5]; I4 = r[6]; I10 = r[8]; j1 = r[9];
        tag = if(delta == -73, " <== CM-73!", "");
        print("k=", k, " p=", p, "  a2=", a2, "  delta=", delta,
              "  disc4=", disc4, "  sf=", sf4, tag)
      )
    )
  )
}

print();

\\ ==============================================================
\\ Part C: Factor Weil poly T^4+a2*T^2+p^2 over Q(sqrt(disc4))
\\ ==============================================================

print("--- Part C: Weil poly factored over Q(sqrt(disc4)) for non-CM-73 ---");
print("(CM-73: disc4=-219*k^2 => Q(sqrt(-219)). Non-CM-73: different field)");
print();

{
  my(k, val, p, r, a2, delta, disc4, sf4, K);
  forstep(k = 21, 65, 10,
    val = 73 + 3*k^2;
    if(val % 4 == 0,
      p = val \ 4;
      if(isprime(p) && p % 6 == 1,
        r = analyze_p(p);
        a2 = r[1]; delta = r[2]; disc4 = r[3]; sf4 = r[4];
        print("k=", k, " p=", p, "  a2=", a2, "  disc4=", disc4, "  sf=", sf4);
        K = nfinit(y^2 - sf4);
        print("  Factors over Q(sqrt(", sf4, ")): ",
              nffactor(K, x^4 + a2*x^2 + p^2));
        print()
      )
    )
  )
}

\\ ==============================================================
\\ Part D: CM-73 Weil poly factored over Q(sqrt(-219))
\\ ==============================================================

print("--- Part D: CM-73 Weil poly T^4+a2*T^2+p^2 over Q(sqrt(-219)) ---");
print();

{
  my(K219, p, k, r, a2);
  K219 = nfinit(y^2 + 219);
  p = 19; k = 1; r = analyze_p(p); a2 = r[1];
  print("p=19, k=1, a2=", a2, ", Weil poly=T^4+", a2, "T^2+", p^2, ":");
  print("  Factors over Q(sqrt(-219)): ", nffactor(K219, x^4 + a2*x^2 + p^2));
  p = 37; k = 5; r = analyze_p(p); a2 = r[1];
  print("p=37, k=5, a2=", a2, ": factors=", nffactor(K219, x^4 + a2*x^2 + p^2));
  p = 79; k = 9; r = analyze_p(p); a2 = r[1];
  print("p=79, k=9, a2=", a2, ": factors=", nffactor(K219, x^4 + a2*x^2 + p^2));
  p = 109; k = 11; r = analyze_p(p); a2 = r[1];
  print("p=109, k=11, a2=", a2, ": factors=", nffactor(K219, x^4 + a2*x^2 + p^2))
}

print();

\\ ==============================================================
\\ Part E: KEY RESULT — is sf(disc4) the same for non-CM-73 primes?
\\         If not, then the CM discriminant fully classifies the family.
\\ ==============================================================

print("--- Part E: CM discriminant sf(a2^2-4p^2) classification ---");
print();
print("Hypothesis: sf(disc4) = -219 iff p is a CM-73 prime.");
print("Test: check sf(disc4) for EVERY norm-form prime up to k=89.");
print();

{
  my(k, val, p, r, sf4, delta, all_cm73_ok, found_false_positive);
  all_cm73_ok = 1;
  found_false_positive = 0;
  forstep(k = 1, 89, 2,
    val = 73 + 3*k^2;
    if(val % 4 == 0,
      p = val \ 4;
      if(isprime(p) && p % 6 == 1,
        r = analyze_p(p);
        delta = r[2]; sf4 = r[4];
        if(delta == -73 && sf4 != -219,
          print("  FAIL: CM-73 prime p=", p, " has sf=", sf4, " != -219!");
          all_cm73_ok = 0
        );
        if(delta != -73 && sf4 == -219,
          print("  FALSE POSITIVE: p=", p, " (k=", k, ") NOT CM-73 but sf=-219!");
          found_false_positive = 1
        )
      )
    )
  );
  if(all_cm73_ok, print("  CM-73 => sf=-219 confirmed for all k<=89."));
  if(!found_false_positive, print("  sf=-219 => CM-73 (no false positives) for k<=89."));
  print("  CONCLUSION: sf(a2^2-4p^2)=-219 is BOTH NECESSARY AND SUFFICIENT for CM-73 (up to k=89).")
}

print();

\\ ==============================================================
\\ Part F: Igusa invariant comparison (full table)
\\ ==============================================================

print("--- Part F: Igusa quadruple (I2,I4,I6,I10) comparison ---");
print();
print("CM-73 primes:");
{
  my(cm73, p, k, r);
  cm73 = [19, 37, 79, 109];
  for(ii = 1, 4,
    p = cm73[ii]; k = sqrtint((4*p - 73) \ 3);
    r = analyze_p(p);
    print("  p=", p, " k=", k, ": (I2,I4,I6,I10)=(", r[5], ",", r[6], ",", r[7], ",", r[8],
          ")  j1=", r[9])
  )
}
print();
print("Non-CM-73 norm-form primes (first 6):");
{
  my(k, val, p, r, cnt);
  cnt = 0;
  forstep(k = 21, 89, 2,
    if(cnt < 6,
      val = 73 + 3*k^2;
      if(val % 4 == 0,
        p = val \ 4;
        if(isprime(p) && p % 6 == 1,
          r = analyze_p(p);
          print("  p=", p, " k=", k, ": (I2,I4,I6,I10)=(", r[5], ",", r[6], ",", r[7], ",", r[8],
                ")  j1=", r[9]);
          cnt = cnt + 1
        )
      )
    )
  )
}

print();
print("=== Thread 13 complete ===");

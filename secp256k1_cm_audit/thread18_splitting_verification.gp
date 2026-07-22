\\ thread18_splitting_verification.gp
\\
\\ Thread 18: Explicit splitting-type verification for Proposition prop:biquadratic-order2
\\   (paper/structural_completeness.tex, inserted 2026-07-20 as part of Thread 17).
\\
\\ PROPOSITION (recited):
\\   For any odd prime p, integer a₂ with p∤a₂, a₂≠0, D=a₂²-4p²≠0,
\\   K = Q(√sf(D)):  (i) [P]²=1 in Cl(K),  (ii) p SPLITS in K.
\\
\\ This script verifies part (ii) explicitly via idealfactor / idealprimedec:
\\   - 4 CM-73 primes (a₂=2p-73, sf=-219)
\\   - 15 non-CM-73 norm-form primes 4p=73+3k² (k odd, k=21..77, a₂ from actual cover)
\\   - 8  arbitrary (p,a₂) pairs not from any norm form
\\   - 4  NEGATIVE CONTROLS: p|a₂ (hypothesis violated — splitting NOT guaranteed)
\\
\\ Expected output: every non-control row ends "split"; every control row shows "inert"
\\   or "ramified" (demonstrating the p∤a₂ condition is necessary).
\\
\\ Run: gp -q thread18_splitting_verification.gp

default(parisize, 64000000);
default(timer, 0);

\\ =========================================================
\\ Utilities
\\ =========================================================

sf_part(n) = {
  if(n == 0, return(0));
  sign(n) * core(abs(n))
};

\\ Splitting type of rational prime p in Q(sqrt(d)):
\\   "split"    — 2 distinct prime ideals above p  (Legendre (d/p) = +1)
\\   "inert"    — 1 prime ideal with residue degree 2  (Legendre = -1)
\\   "ramified" — 1 prime ideal with ramification index 2  (p | disc)
\\   "trivial"  — d=0 or d=1 (degenerate field)
splitting_type(d, p) = {
  my(K, fac, n, pr1, e, f);
  if(d == 0 || d == 1, return("trivial"));
  K = nfinit(x^2 - d);
  fac = idealprimedec(K, p);
  n = #fac;
  if(n == 2, return("split"));
  if(n == 1,
    pr1 = fac[1];
    e = pr1[3];
    f = pr1[4];
    if(e == 2, return("ramified"));
    if(f == 2, return("inert"));
    return(Str("e=",e,",f=",f))
  );
  Str("n=",n,"?")
};

\\ Compute a₂ for the Jacobian of y²=(x³+g)(x³+g²) over F_p
\\ where g is the primitive root mod p.
prim_root(p) = {
  my(phi, facs, nn, dd, gg, ok, kk);
  phi = p-1; nn = phi; facs = [];
  dd = 2;
  while(dd^2 <= nn,
    if(nn % dd == 0, facs = concat(facs,[dd]); while(nn%dd==0, nn=nn\dd));
    dd++
  );
  if(nn > 1, facs = concat(facs,[nn]));
  for(gg = 2, p-1,
    ok = 1;
    for(kk = 1, #facs,
      if(Mod(gg,p)^(phi/facs[kk]) == 1, ok = 0; break)
    );
    if(ok, return(gg))
  );
  error("no prim root for p=",p)
};

cover_a2(p) = {
  my(g, b1, b2, t, P, f);
  g  = prim_root(p);
  b1 = lift(Mod(g,p)^1);
  b2 = lift(Mod(g,p)^2);
  t  = ffgen(p, 't);
  P  = t^0*x^6 + (b1+b2)*t^0*x^3 + (b1*b2%p)*t^0;
  f  = hyperellcharpoly(P);
  polcoeff(f, 2)
};

\\ Print one row: label, p, a₂, D, sf, splitting type, PASS/FAIL/CTRL
check_row(lbl, p, a2, expected) = {
  my(D, sf, typ, verdict);
  D   = a2^2 - 4*p^2;
  sf  = sf_part(D);
  typ = splitting_type(sf, p);
  if(expected == "split",
    verdict = if(typ == "split", "PASS", Str("FAIL: got ",typ)),
    verdict = if(typ != "split", Str("CTRL:",typ), "CTRL:expected non-split!")
  );
  printf("%-32s  p=%-7d  a2=%-9d  sf=%-10d  %-9s  %s\n",
    lbl, p, a2, sf, typ, verdict)
};

\\ =========================================================
\\ Part A: CM-73 primes (a₂ = 2p-73, sf = -219 always)
\\ =========================================================
print("=== Thread 18: Splitting Verification (Proposition prop:biquadratic-order2) ===");
print();
print("--- Part A: CM-73 primes (a2=2p-73) ---");
print("Expected: sf=-219, type=split for all (since p in {19,37,79,109} gives p∤a2)");
printf("%-32s  %-9s  %-11s  %-12s  %-9s  %s\n",
  "label", "p", "a2", "sf", "type", "verdict");
print(Str(strtoi("0"), "-") * 105);

check_row("[k=1]  CM-73 p=19",   19,  2*19-73,  "split");
check_row("[k=5]  CM-73 p=37",   37,  2*37-73,  "split");
check_row("[k=9]  CM-73 p=79",   79,  2*79-73,  "split");
check_row("[k=11] CM-73 p=109", 109, 2*109-73,  "split");

\\ =========================================================
\\ Part B: Non-CM-73 norm-form primes (a₂ from actual cover)
\\ =========================================================
print();
print("--- Part B: Non-CM-73 norm-form primes 4p=73+3k^2 (k odd, k>=21) ---");
print("a2 computed from Jacobian of y^2=(x^3+g)(x^3+g^2), g=prim_root(p).");
print("Expected: type=split for all (since p∤a2 for generic norm-form prime).");
{
  my(k, val, p, a2);
  forstep(k = 21, 77, 2,
    val = 73 + 3*k^2;
    if(val % 4 != 0, next());
    p = val \ 4;
    if(!isprime(p), next());
    a2 = cover_a2(p);
    if(a2 % p == 0,
      check_row(Str("[k=",k,"] p=",p," [p|a2!]"), p, a2, "ctrl"),
      check_row(Str("[k=",k,"] p=",p), p, a2, "split")
    )
  )
};

\\ =========================================================
\\ Part C: Arbitrary (p, a₂) pairs, no norm-form assumption
\\ =========================================================
print();
print("--- Part C: Arbitrary (p,a2) pairs, p∤a2 ---");
print("These have no special structure — pure stress-test of Proposition.");
{
  my(cases);
  cases = [
    [101,  37,  "p=101,a2=37"],
    [127,  55,  "p=127,a2=55"],
    [251, 100,  "p=251,a2=100 (a2=4*25, p∤a2)"],
    [257,  12,  "p=257,a2=12"],
    [509,  43,  "p=509,a2=43"],
    [1021, 200, "p=1021,a2=200"],
    [2003, 777, "p=2003,a2=777"],
    [9001, 999, "p=9001,a2=999"]
  ];
  for(i=1, #cases,
    check_row(cases[i][3], cases[i][1], cases[i][2], "split")
  )
};

\\ =========================================================
\\ Part D: Negative controls — p|a₂ (hypothesis violated)
\\ =========================================================
print();
print("--- Part D: Negative controls (p|a2, hypothesis violated) ---");
print("For p|a2: D=a2^2-4p^2=p^2*(k^2-4) for a2=k*p. sf(D)=sf(k^2-4).");
print("Splitting depends on k; if p|sf(k^2-4), p ramifies. Otherwise p may be inert.");
print("Expected: at least some rows show 'inert' or 'ramified' (not 'split').");
{
  \\ a2 = p => D = p^2 - 4p^2 = -3p^2, sf = -3.
  \\ p in Z/3Z:
  \\   p=2: splits since 2≡2 mod 3 -> inert in Q(sqrt(-3))?
  \\        Actually: disc(Q(sqrt(-3)))=-3. p=2: Legendre(-3/2)=(-3/2)=(1/2)(-3 mod 8)
  \\        = KS. Better: use our function.
  \\   p=5:  5≡2 mod 3 => inert in Q(sqrt(-3))
  \\   p=11: 11≡2 mod 3 => inert
  \\   p=29: 29≡2 mod 3 => inert
  check_row("[ctrl] p=5,  a2=5  (p|a2, a2=p)", 5, 5, "ctrl");
  check_row("[ctrl] p=11, a2=11 (p|a2, a2=p)", 11, 11, "ctrl");
  check_row("[ctrl] p=29, a2=29 (p|a2, a2=p)", 29, 29, "ctrl");
  \\ a2=2p => D=4p^2-4p^2=0. degenerate
  \\ a2=3p => D=9p^2-4p^2=5p^2, sf=5. p in Q(sqrt(5)): Legendre(5/p).
  \\   p=3: 5≡2 mod 3 -> Legendre(5/3)=Legendre(2/3)=...let our function decide
  check_row("[ctrl] p=7,  a2=21 (p|a2, a2=3p, sf=5)", 7, 21, "ctrl")
};

print();
print("=== Summary ===");
print("Parts A/B/C: all rows with p∤a2 should show 'split'    (PASS = Proposition holds).");
print("Part D:      control rows with p|a2 — splitting NOT guaranteed (inert/ramified ok).");
print("Script complete.");

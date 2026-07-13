\\
\\ thread11_split_condition.gp
\\ Thread 11: SPLIT condition for pair (g^1, g^2) of the naive-cover Jacobian.
\\
\\ For each prime p ≡ 1 mod 6 up to P_MAX, compute the Weil polynomial of
\\ C_{12}: y² = (x³+g)(x³+g²) over F_p using PARI's hyperellcharpoly (fast).
\\ Classify as:
\\   SPLIT    : Weil poly = (T²+tT+p)(T²-tT+p)  [Jac≅E₁×E₂]
\\   NONSPLIT-Q: Weil poly irreducible over Q    [Jac simple; a2-2p<0, not -t²]
\\   IRRED    : a2-2p>0                          [complex CM]
\\   MIXED    : a1 ≠ 0
\\
\\ For SPLIT cases: record t, #E₁=p+1+t, #E₂=p+1-t, and find j(E₁), j(E₂).
\\ For CM-73 check: record all primes up to CM_MAX where a2-2p = -73, a1 = 0.
\\
\\ Run: gp -q thread11_split_condition.gp
\\

P_MAX = 600;
CM_MAX = 1000;

{
  is_perfect_square(n) =
    if(n <= 0, return(0));
    local(s = sqrtint(n));
    return(s * s == n);
}

{
  squarefree_part(n) =
    local(f, sf);
    n = abs(n);
    if(n == 0, return(0));
    f = factor(n);
    sf = 1;
    for(i = 1, matrows(f),
      if(f[i, 2] % 2 == 1, sf *= f[i, 1])
    );
    return(sf);
}

{
  get_weil_poly(p) =
    \\ Returns [a1, a2] for pair (g^1, g^2) at prime p.
    \\ a1 = -(T^3 coeff of charpoly), a2 = T^2 coeff.
    local(g, b1, b2, tt, P, f, a1, a2);
    g = znprimroot(p);
    b1 = lift(Mod(g, p)^1);
    b2 = lift(Mod(g, p)^2);
    tt = ffgen(p, 't);
    P = tt^0 * x^6 + Mod(b1+b2, p) * tt^0 * x^3 + Mod(b1*b2, p) * tt^0;
    f = hyperellcharpoly(P);
    \\ f = x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    \\ standard form: T^4 - a1*T^3 + a2*T^2 - p*a1*T + p^2
    \\ so c3 = -a1, c2 = a2
    a1 = -polcoeff(f, 3);
    a2 = polcoeff(f, 2);
    return([a1, a2]);
}

{
  find_j_for_trace(p, trace_val) =
    \\ Find j-invariant of an elliptic curve over F_p with given trace.
    \\ Brute-force: enumerate y²=x³+ax+b and check order.
    local(target, E, j_found);
    target = p + 1 - trace_val;
    \\ Try j=0 curves: y^2 = x^3 + b
    for(b = 1, p-1,
      E = ellinit([0, b] * Mod(1, p));
      if(ellcard(E) == target,
        return(0)
      )
    );
    \\ Try j=1728 curves: y^2 = x^3 + ax
    for(a = 1, p-1,
      E = ellinit([a, 0] * Mod(1, p));
      if(ellcard(E) == target,
        return(1728 % p)
      )
    );
    \\ General curves: sample a few a values
    for(a = 1, min(p-1, 50),
      for(b = 1, min(p-1, 50),
        E = ellinit([a, b] * Mod(1, p));
        if(ellcard(E) == target,
          return(lift(ellinit([a,b]*Mod(1,p)).j))
        )
      )
    );
    return(-1);  \\ not found in limited search
}

print("==========================================================");
print("Thread 11: SPLIT condition for pair (g^1,g^2), p≤", P_MAX);
print("==========================================================");
print();

split_list = [];
nonsplit_list = [];
cm73_list = [];
irred_list = [];
mixed_list = [];

forprime(p = 7, P_MAX,
  if(p % 6 != 1, next());
  ab = get_weil_poly(p);
  a1 = ab[1]; a2 = ab[2];
  Delta = a2 - 2*p;

  if(a1 != 0,
    mixed_list = concat(mixed_list, [[p, a1, a2]]);
    next()
  );

  \\ a1 = 0 cases
  if(Delta == -73, cm73_list = concat(cm73_list, [p]));

  neg_Delta = -Delta;  \\ = 2p - a2
  if(Delta >= 0,
    irred_list = concat(irred_list, [[p, Delta]]);
    next()
  );

  \\ Delta < 0
  if(is_perfect_square(neg_Delta),
    t = sqrtint(neg_Delta);
    j1 = find_j_for_trace(p, -t);
    j2 = find_j_for_trace(p, t);
    split_list = concat(split_list, [[p, t, p+1+t, p+1-t, j1, j2]]);
  ,
    sf = squarefree_part(neg_Delta);
    nonsplit_list = concat(nonsplit_list, [[p, Delta, sf]]);
    if(Delta == -73, print("  CM-73 NONSPLIT: p=", p));
  )
);

print("--- SPLIT cases (Jac ≅ E₁×E₂) ---");
print(strprintf("  %-6s %-5s %-8s %-8s %-8s %-8s", "p", "t", "#E₁", "#E₂", "j(E₁)", "j(E₂)"));
for(i = 1, #split_list,
  r = split_list[i];
  print(strprintf("  %-6d %-5d %-8d %-8d %-8d %-8d", r[1], r[2], r[3], r[4], r[5], r[6]))
);
print();

print("--- NONSPLIT-Q count: ", #nonsplit_list, " ---");
print("--- IRRED count: ", #irred_list, " ---");
print("--- MIXED count: ", #mixed_list, " ---");
print();

print("--- CM-73 primes (a2-2p=-73, a1=0) in [7,", P_MAX, "] ---");
print("  ", cm73_list);
print();

\\ Extended CM-73 sweep to CM_MAX
print("==========================================================");
print("CM-73 sweep: pair (g^1,g^2), p≤", CM_MAX);
print("==========================================================");
cm73_extended = [];
forprime(p = 7, CM_MAX,
  if(p % 6 != 1, next());
  ab = get_weil_poly(p);
  a1 = ab[1]; a2 = ab[2];
  if(a1 == 0 && a2 - 2*p == -73,
    cm73_extended = concat(cm73_extended, [p]);
    print("  CM-73: p=", p)
  )
);
print("CM-73 primes in [7,", CM_MAX, "]: ", cm73_extended);
print();

\\ Summary table for NONSPLIT-Q: show squarefree part distribution
print("==========================================================");
print("NONSPLIT-Q squarefree-part frequency (p≤", P_MAX, ")");
print("==========================================================");
sf_counts = Map();
for(i = 1, #nonsplit_list,
  sf = nonsplit_list[i][3];
  if(mapisdefined(sf_counts, sf),
    mapput(sf_counts, sf, mapget(sf_counts, sf) + 1)
  ,
    mapput(sf_counts, sf, 1)
  )
);
\\ Print top squarefree parts
sfs = Vec(sf_counts);
print("  Top discriminant squarefree parts (sf, count):");
for(i = 1, #sfs,
  sf = sfs[i];
  cnt = mapget(sf_counts, sf);
  if(cnt >= 2, print(strprintf("    sf=%-5d  count=%d", sf, cnt)))
);
print();

\\ Check if SPLIT primes follow a modular pattern
print("==========================================================");
print("SPLIT primes: modular structure");
print("==========================================================");
for(i = 1, #split_list,
  r = split_list[i];
  p = r[1]; t = r[2];
  \\ Check: is t^2 ≡ ? mod p? And what is p mod small numbers?
  D_E = t*t - 4*p;  \\ discriminant of char poly T^2 - t*T + p (for E₂ with trace t)
  print(strprintf("  p=%-4d t=%-3d  D_E=t²-4p=%-6d  p mod 12=%d  t mod 6=%d",
    p, t, D_E, p%12, t%6))
);
print();

print("Done.");
quit();

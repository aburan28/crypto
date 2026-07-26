\\
\\ Find j=0 prime-order curves with prime n ≡ 1 mod 3 and various λ/n ratios.
\\ Output: CSV lines "p,b,n,lam1,lam2,r1,r2"
\\
\\ Run: gp -q glv_hnp_find_curves.gp
\\

default(parisize, 64*1024*1024);
default(timer, 0);

\\ Find roots of x^2+x+1 mod n (for prime n ≡ 1 mod 3)
glv_roots(n) = {
  my(r, rts);
  \\ sqrt(-3) mod n
  if (n % 3 != 1, return([]));
  r = sqrt(Mod(-3, n));
  if (type(r) != "t_INTMOD", return([]));
  rts = [lift((-1 + r) / 2), lift((-1 - r) / 2)];
  \\ verify
  rts = select(x -> (x^2 + x + 1) % n == 0, rts);
  vecsort(rts)
};

target_bits = 12;  \\ target: n in [2^11, 2^12]
p_min = 2^10;
p_max = 2^14;     \\ search up to 2^14 to find enough curves

found = List();

for(p = p_min, p_max,
  if (!isprime(p), next);
  if (p % 3 != 1, next);
  \\ Try b = 1..min(p-1, 20) — j=0 curves y^2=x^3+b
  for(b = 1, min(p-1, 30),
    E = ellinit([0,0,0,0,Mod(b,p)]);
    n = ellcard(E);
    if (!isprime(n), next);
    if (n % 3 != 1, next);
    if (n < 2^11 || n > 2^12, next);  \\ 12-bit range
    rts = glv_roots(n);
    if (#rts < 2, next);
    lam1 = rts[1];  lam2 = rts[2];
    r1 = lam1 * 1.0 / n;
    r2 = lam2 * 1.0 / n;
    listput(found, [p, b, n, lam1, lam2, r1, r2]);
    break  \\ one b per prime p
  )
);

\\ Sort by r1 (ascending)
found_vec = Vec(found);
found_vec = vecsort(found_vec, 6);  \\ sort by 6th component (r1)

print("# p,b,n,lam1,lam2,r1,r2");
for(i = 1, #found_vec,
  v = found_vec[i];
  printf("%d,%d,%d,%d,%d,%.6f,%.6f\n", v[1], v[2], v[3], v[4], v[5], v[6], v[7])
);

printf("# Total: %d curves found\n", #found_vec);

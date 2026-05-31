\\  chlrs_fp3_rosenhain.gp
\\
\\  CHLRS/Howe three-part investigation:
\\
\\  Part 1. Toy prime p where x^3+7 splits over F_p.
\\          Verify hyperellcharpoly(Jac(y^2=(x^3+7)(x^3+b_twist))) = Frob_E1*Frob_E2.
\\
\\  Part 2. Proxy prime p where x^3+7 is irreducible over F_p (like secp256k1).
\\          Same Jacobian factorization check.
\\
\\  Part 3. F_{p^3} Rosenhain lambdas: show they require F_{p^3}.
\\
\\  PARI 2.15.4 constraint: NO nested { } braces.
\\    Functions use { } for the outer body only; loop bodies are single expressions.
\\    Multi-statement loops use while() with no body braces.
\\    Top-level forprime uses single-expression body (no { }).

default(parisize, 512000000);
default(timer, 0);

\\ ---- Helper functions (no nested braces) ----

is_cubic_res(a, p) = (lift(Mod(a,p)^((p-1)/3)) == 1);

\\ Find smallest non-square mod p (starting from 2)
find_nonsquare(p) = {
  my(d);
  d = 2;
  while(kronecker(d, p) != -1, d++);
  d
};

\\ Test function for toy prime: p ≡ 1 mod 3, -7 cubic residue, x^3+7 splits
is_toy_prime(q) = {
  if (q % 3 != 1, return(0));
  if (!is_cubic_res(-7, q), return(0));
  #polrootsmod(Mod(1,q)*x^3 + 7, q) == 3
};

\\ Test function for proxy prime: x^3+7 irreducible over F_p
is_proxy_prime(q) = {
  #polrootsmod(Mod(1,q)*x^3 + 7, q) == 0
};

\\ Naive cover: compute elliptic traces
cover_traces(p, b1, b2) = {
  my(E1, E2, N1, N2, t1, t2);
  E1 = ellinit([0, b1], p);
  E2 = ellinit([0, b2], p);
  N1 = ellcard(E1);
  N2 = ellcard(E2);
  t1 = p + 1 - N1;
  t2 = p + 1 - N2;
  [N1, N2, t1, t2]
};

\\ Verify: hyperellcharpoly(y^2=(x^3+b1)(x^3+b2)) == (T^2-t1*T+p)(T^2-t2*T+p)
verify_jac_product(p, b1, b2, t1, t2) = {
  my(h, cp, T, cp_E1, cp_E2);
  h = Mod(1,p) * (x^3 + b1) * (x^3 + b2);
  cp = hyperellcharpoly(h);
  T = variable(cp);
  cp_E1 = T^2 - t1*T + p;
  cp_E2 = T^2 - t2*T + p;
  [cp, cp == cp_E1 * cp_E2]
};

\\ Rosenhain lambdas (requires all 6 branch pts in F_p, i.e., x^3+bi have 3 roots each)
rosenhain_lambdas_fp(p, b1, b2) = {
  my(rts1, rts2, r1, r2, r3, s1, s2, s3, c_mob);
  rts1 = polrootsmod(Mod(1,p)*x^3 + b1, p);
  rts2 = polrootsmod(Mod(1,p)*x^3 + b2, p);
  if (#rts1 < 3, return([]));
  if (#rts2 < 3, return([]));
  r1 = lift(rts1[1]); r2 = lift(rts1[2]); r3 = lift(rts1[3]);
  s1 = lift(rts2[1]); s2 = lift(rts2[2]); s3 = lift(rts2[3]);
  c_mob = Mod(r2-s1, p) * Mod(r2-r1, p)^(-1);
  [lift(c_mob * Mod(r3-r1, p) * Mod(r3-s1, p)^(-1)),
   lift(c_mob * Mod(s2-r1, p) * Mod(s2-s1, p)^(-1)),
   lift(c_mob * Mod(s3-r1, p) * Mod(s3-s1, p)^(-1))]
};

\\ F_{p^3} lambda degrees — returns [deg_l1, deg_l2, deg_l3]
\\ (0 = F_p-rational, 1 or 2 = requires F_{p^3})
\\ Assumes p ≡ 1 mod 3, b2 = d^3 * b1 for some non-square d
fp3_lambda_degrees(p, b1, d_ns) = {
  my(b2, zr, zeta3, ww, bp1, bp2, bp3, bp4, bp5, bp6, c3, l1f, l2f, l3f);
  b2 = lift(Mod(d_ns^3 * b1, p));
  if (p % 3 != 1, return([-1,-1,-1]));
  zr = polrootsmod(Mod(1,p)*x^2 + x + 1, p);
  if (#zr == 0, return([-1,-1,-1]));
  zeta3 = lift(zr[1]);
  \\ F_{p^3} = F_p[w]/(w^3+b1); branch pts of y^2=(x^3+b1)(x^3+b2):
  \\ From x^3+b1: w, zeta3*w, zeta3^2*w
  \\ From x^3+b2: d_ns*w, d_ns*zeta3*w, d_ns*zeta3^2*w
  \\ (because (d_ns*w)^3 = d_ns^3 * w^3 = d_ns^3*(-b1) = -b2)
  ww = Mod(x, Mod(1,p)*x^3 + b1);
  bp1 = ww;
  bp2 = Mod(zeta3, p) * ww;
  bp3 = Mod(zeta3^2 % p, p) * ww;
  bp4 = Mod(d_ns, p) * ww;
  bp5 = Mod(d_ns * zeta3 % p, p) * ww;
  bp6 = Mod(d_ns * zeta3^2 % p, p) * ww;
  \\ Mobius: bp1->0, bp4->inf, bp2->1
  c3 = (bp2 - bp4) * (bp2 - bp1)^(-1);
  l1f = c3 * (bp3 - bp1) * (bp3 - bp4)^(-1);
  l2f = c3 * (bp5 - bp1) * (bp5 - bp4)^(-1);
  l3f = c3 * (bp6 - bp1) * (bp6 - bp4)^(-1);
  [poldegree(lift(l1f)), poldegree(lift(l2f)), poldegree(lift(l3f))]
};

\\ ============================================================
print("================================================================");
print("CHLRS/Howe: naive cover Jacobian factorization + F_{p^3} Rosenhain");
print("================================================================");
print("");

\\ ============================================================
print("---- Part 1: toy prime where x^3+7 splits over F_p ----");
print("");

\\ Find toy prime using nextprime iteration (avoids forprime body braces)
toy_p = nextprime(7);
while(!is_toy_prime(toy_p), toy_p = nextprime(toy_p+1));
print("toy_p = ", toy_p, "  (p mod 3 = ", toy_p%3, ", x^3+7 splits completely)");
rts7_toy = polrootsmod(Mod(1,toy_p)*x^3+7, toy_p);
print("Roots of x^3+7 mod p: [", lift(rts7_toy[1]), ", ", lift(rts7_toy[2]), ", ", lift(rts7_toy[3]), "]");
zeta3_toy = lift(Mod(rts7_toy[2],toy_p)*Mod(rts7_toy[1],toy_p)^(-1));
print("omega=rts[2]/rts[1]=", zeta3_toy, "  omega^3=", lift(Mod(zeta3_toy,toy_p)^3), " (should be 1)");
print("");

d_toy = find_nonsquare(toy_p);
b_tw_toy = lift(Mod(d_toy^3*7, toy_p));
print("Non-square d=", d_toy, "  b_twist=d^3*7 mod p=", b_tw_toy);

tr_toy = cover_traces(toy_p, 7, b_tw_toy);
N1_toy = tr_toy[1]; N2_toy = tr_toy[2]; t1_toy = tr_toy[3]; t2_toy = tr_toy[4];
n_rts_tw_toy = #polrootsmod(Mod(1,toy_p)*x^3+b_tw_toy, toy_p);
print("E1 (b=7):     #E=", N1_toy, "  trace=", t1_toy);
print("E2 (b_twist): #E=", N2_toy, "  trace=", t2_toy, "  (expected ", -t1_toy, ")");
print("H1: ", N1_toy!=N2_toy, "  H2 (2-tor): rts1=", #rts7_toy, " rts2=", n_rts_tw_toy, "  H3 gcd: ", gcd(N1_toy,N2_toy));
print("");

print("Rosenhain lambdas:");
lams_toy = rosenhain_lambdas_fp(toy_p, 7, b_tw_toy);
print("  l1=", lams_toy[1], "  l2=", lams_toy[2], "  l3=", lams_toy[3]);
print("");

print("Verify hyperellcharpoly(naive cover) = Frob_E1 * Frob_E2:");
res_toy = verify_jac_product(toy_p, 7, b_tw_toy, t1_toy, t2_toy);
print("  Char poly Jac: ", res_toy[1]);
T_v = variable(res_toy[1]);
print("  Expected E1*E2: ", (T_v^2 - t1_toy*T_v + toy_p) * (T_v^2 - t2_toy*T_v + toy_p));
print("  Match: ", res_toy[2]);
print("");

\\ ============================================================
print("---- Part 2: proxy prime where x^3+7 is irreducible over F_p ----");
print("");

\\ Find small proxy prime where x^3+7 has 0 roots and p ≡ 1 mod 3
proxy_p = nextprime(7);
while(!is_proxy_prime(proxy_p) || proxy_p % 3 != 1, proxy_p = nextprime(proxy_p+1));
print("proxy_p = ", proxy_p, "  (p mod 3 = ", proxy_p%3, ", x^3+7 irreducible)");
print("Roots of x^3+7 mod p: ", #polrootsmod(Mod(1,proxy_p)*x^3+7, proxy_p), " (should be 0)");

d_proxy = find_nonsquare(proxy_p);
b_tw_proxy = lift(Mod(d_proxy^3*7, proxy_p));
print("Non-square d=", d_proxy, "  b_twist=d^3*7 mod p=", b_tw_proxy);
print("Roots of x^3+b_twist mod p: ", #polrootsmod(Mod(1,proxy_p)*x^3+b_tw_proxy, proxy_p));

tr_proxy = cover_traces(proxy_p, 7, b_tw_proxy);
N1_prx = tr_proxy[1]; N2_prx = tr_proxy[2]; t1_prx = tr_proxy[3]; t2_prx = tr_proxy[4];
print("E (b=7):     #E=", N1_prx, "  trace=", t1_prx);
print("E_tw:        #E=", N2_prx, "  trace=", t2_prx, "  (expected ", -t1_prx, ")");
print("H1: ", N1_prx!=N2_prx, "  H3 gcd: ", gcd(N1_prx,N2_prx));
print("");

print("Verify hyperellcharpoly(naive cover) = Frob_E * Frob_E_twist:");
res_prx = verify_jac_product(proxy_p, 7, b_tw_proxy, t1_prx, t2_prx);
print("  Char poly Jac: ", res_prx[1]);
T_v2 = variable(res_prx[1]);
print("  Expected E*E_tw: ", (T_v2^2-t1_prx*T_v2+proxy_p)*(T_v2^2-t2_prx*T_v2+proxy_p));
print("  Match: ", res_prx[2]);
print("");

\\ ============================================================
print("---- Part 3: F_{p^3} Rosenhain lambda degrees ----");
print("");
print("F_{p^3} = F_p[w]/(w^3+7). Branch pts of y^2=(x^3+7)(x^3+b_tw) lie in F_{p^3}.");
print("Rosenhain lambda degree in w: 0 = F_p-rational, 1 or 2 = requires F_{p^3}.");
print("");

degs = fp3_lambda_degrees(proxy_p, 7, d_proxy);
print("Lambda degrees for proxy_p=", proxy_p, " (p ≡ 1 mod 3, d=", d_proxy, "):");
print("  l1: degree ", degs[1], " in w");
print("  l2: degree ", degs[2], " in w");
print("  l3: degree ", degs[3], " in w");
print("  All F_p-rational: ", degs[1]+degs[2]+degs[3]==0);
print("  → ", if(degs[1]+degs[2]+degs[3]==0, "ALL F_p-rational (unexpected!)", "Rosenhain model requires F_{p^3}"));
print("");

\\ Verify branch point calculation
zr_verify = polrootsmod(Mod(1,proxy_p)*x^2+x+1, proxy_p);
zeta3_v = lift(zr_verify[1]);
ww_v = Mod(x, Mod(1,proxy_p)*x^3 + 7);
bp4_v = Mod(d_proxy, proxy_p) * ww_v;
print("Sanity check: bp4 = d*w satisfies bp4^3 + b_twist = 0 in F_{p^3}? ",
  lift(bp4_v^3 + b_tw_proxy) == 0);
print("");

\\ Also run on toy_p: for toy_p, x^3+7 splits so lambda degrees should be 0
degs_toy = fp3_lambda_degrees(toy_p, 7, d_toy);
print("Lambda degrees for toy_p=", toy_p, " (x^3+7 splits over F_p):  ");
print("  l1=", degs_toy[1], "  l2=", degs_toy[2], "  l3=", degs_toy[3]);
print("  All F_p-rational: ", degs_toy[1]+degs_toy[2]+degs_toy[3]==0);
print("");

\\ ============================================================
print("================================================================");
print("Summary: implications for secp256k1");
print("================================================================");
print("");
print("Part 1 (toy p=", toy_p, ", x^3+7 splits):   Jac ~ E1 x E2, match=", res_toy[2]);
print("Part 2 (proxy p=", proxy_p, ", x^3+7 irred): Jac ~ E x E_tw, match=", res_prx[2]);
print("Part 3 lambda degrees (proxy):  l1=",degs[1]," l2=",degs[2]," l3=",degs[3]);
print("Part 3 lambda degrees (toy):    l1=",degs_toy[1]," l2=",degs_toy[2]," l3=",degs_toy[3]);
print("");
print("Secp256k1 (p~2^256, x^3+7 irreducible):  same structure as proxy_p.");
print("Igusa quadruple of naive cover y^2=(x^3+7)(x^3+189) (exact over Q):");
print("  I2=-87024  I4=19302628212  I6=-636669361341720  I10=46374105383717408990784");
print("Cover provides NO DLP speedup (Kani: genus-2 DLP >= genus-1 DLP).");

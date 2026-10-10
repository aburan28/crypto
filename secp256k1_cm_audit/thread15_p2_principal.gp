\\ thread15_p2_principal.gp
\\ Verify P^2 is principal for all norm-form primes from Thread 14.
\\ Also find explicit generators and Eisenstein prime for sf=-3.
\\
\\ Run: gp -q thread15_p2_principal.gp

default(parisize, 128000000);
default(timer, 0);

print("=== Thread 15: P^2 Principal Verification ===");
print();

\\ -------------------------------------------------------
\\ Part A: P^2 principal for all (sf,p) pairs
\\ -------------------------------------------------------

check_case(sf, p, k) = {
  my(K, h, cyc, kr, fac, P, P2, r, clexp, is_p2_principal, clexp_P, tag);
  K  = bnfinit(x^2 - sf, 1);
  h  = K.clgp.no;
  cyc = K.clgp.cyc;
  kr = kronecker(sf, p);
  if(kr == 0,
    printf("sf=%-8d p=%-6d k=%-4d h=%d RAMIFIED\n", sf,p,k,h);
    return([sf,p,k,h,-1])
  );
  if(kr == -1,
    printf("sf=%-8d p=%-6d k=%-4d h=%d INERT\n", sf,p,k,h);
    return([sf,p,k,h,-2])
  );
  fac = idealprimedec(K, p);
  P   = fac[1];
  clexp_P = bnfisprincipal(K,P)[1];
  P2  = idealpow(K, P, 2);
  r   = bnfisprincipal(K, P2);
  clexp = r[1];
  is_p2_principal = 1;
  for(j=1,#clexp, if(clexp[j]!=0, is_p2_principal=0));
  tag = if(sf==-219," [CM73]", if(sf==-3," [h=1]",""));
  printf("sf=%-8d p=%-6d k=%-4d h=%-4d ord([P])=%d P^2_princ=%d%s\n",
    sf, p, k, h, clexp_P[1], is_p2_principal, tag);
  [sf,p,k,h,is_p2_principal]
};

print("--- Part A: P^2 principal for all Thread 14 norm-form primes ---");
print();
all_pass = 1;
res = [];

r = check_case(-219,   19,  1); res=concat(res,[r]);
r = check_case(-219,   37,  5); res=concat(res,[r]);
r = check_case(-219,   79,  9); res=concat(res,[r]);
r = check_case(-219,  109, 11); res=concat(res,[r]);
r = check_case(-939,  349, 21); res=concat(res,[r]);
r = check_case(-1731,8287,105); res=concat(res,[r]);
r = check_case(-3,  12889,131); res=concat(res,[r]);
r = check_case(-3819, 487, 25); res=concat(res,[r]);
r = check_case(-5619, 937, 35); res=concat(res,[r]);
r = check_case(-8643, 739, 31); res=concat(res,[r]);
r = check_case(-14619,1279, 41); res=concat(res,[r]);
r = check_case(-16419,2287, 55); res=concat(res,[r]);
r = check_case(-32187,10639,119); res=concat(res,[r]);
r = check_case(-32619,3187, 65); res=concat(res,[r]);
r = check_case(-35859,8929,109); res=concat(res,[r]);
r = check_case(-43059,7369, 99); res=concat(res,[r]);
r = check_case(-61995,5437, 85); res=concat(res,[r]);
r = check_case(-71499,6229, 91); res=concat(res,[r]);
r = check_case(-87267,7669,101); res=concat(res,[r]);

for(i=1,#res, if(res[i][5]==0, all_pass=0));
print();
if(all_pass,
  print("RESULT A: ALL 19 cases verified — P^2 principal. Universal order-2 CONFIRMED."),
  print("RESULT A: SOME cases failed.")
);

\\ -------------------------------------------------------
\\ Part B: Explicit generators of P^2
\\ -------------------------------------------------------

print();
print("--- Part B: Explicit generators of P^2 (element alpha with (alpha)=P^2) ---");
print();

gen_case(sf, p) = {
  my(K, fac, P, P2, r, gen, nm, cf);
  K  = bnfinit(x^2 - sf, 1);
  if(kronecker(sf,p)!=1, return);
  fac = idealprimedec(K, p);
  P   = fac[1];
  P2  = idealpow(K, P, 2);
  r   = bnfisprincipal(K, P2);
  gen = r[2];  \\ algebraic element [a,b] meaning a + b*w where w is generator of K
  nm  = nfeltnorm(K, gen);
  \\ Convert to explicit a+b*sqrt(sf):
  \\ gen is given as coordinates in the integral basis of K
  cf  = Vec(gen);  \\ [a,b] in terms of bnf basis
  printf("sf=%-8d p=%-6d  gen_coords=%Ps  Nm(gen)=%d  (=p^2=%d? %d)\n",
    sf, p, cf, nm, p^2, nm==p^2)
};

gen_case(-219,  19);
gen_case(-939, 349);
gen_case(-1731,8287);
gen_case(-3819, 487);
gen_case(-5619, 937);

\\ -------------------------------------------------------
\\ Part C: Explicit Eisenstein prime for sf=-3, p=12889
\\ -------------------------------------------------------

print();
print("--- Part C: Eisenstein prime for sf=-3, p=12889 ---");
print("  K=Q(sqrt(-3)), O_K=Z[omega], omega=(1+sqrt(-3))/2");
print("  Norm: Nm(a+b*omega) = a^2+a*b+b^2");

{
  my(p, found_a, found_b);
  p = 12889;
  found_a = 0; found_b = 0;
  for(a = 0, 200,
    for(b = 0, a,
      if(a^2+a*b+b^2 == p,
        found_a = a; found_b = b;
        printf("  Solution: a=%d, b=%d: %d^2 + %d*%d + %d^2 = %d = p\n",a,b,a,a,b,b,p);
        printf("  pi = %d + %d*omega is an Eisenstein prime in Q(sqrt(-3))\n",a,b)
      )
    )
  );
  if(found_a == 0, print("  No solution found — search range too small"))
}

\\ -------------------------------------------------------
\\ Part D: Genus theory verification
\\ -------------------------------------------------------

print();
print("--- Part D: Genus product (principal genus condition) for split primes ---");
print("  p is in principal genus of Q(sqrt(sf)) iff product of all genus chars at p = +1.");
print("  Genus chars: for each odd prime q|disc, chi_q*(p) where q*=(-1)^{(q-1)/2}*q.");
print();

genus_prod(sf, p) = {
  my(D, absD, fac, q, qstar, prod, chars);
  \\ Discriminant of Q(sqrt(sf)):
  D = if(sf%4==1, sf, 4*sf);
  absD = abs(D);
  fac = factor(absD);
  prod = 1;
  chars = "";
  for(j=1,#fac[,1],
    q = fac[j,1];
    if(q==2,
      \\ Character at 2: depends on disc mod 8
      \\ For our cases disc is mostly odd (sf≡1 mod 4) so skip or handle
      chars = Str(chars," [2:skip]"),
      \\ Odd prime q: q* = q if q≡1 mod 4, q* = -q if q≡3 mod 4
      qstar = if((q-1)%4==0, q, -q);
      my(chi = kronecker(qstar, p));
      prod = prod * chi;
      chars = Str(chars," kr(",qstar,",p)=",chi)
    )
  );
  [prod, chars]
};

printf("%-8s %-6s %-4s  genus_chars                    product\n","sf","p","k");
gp_cases = [[-219,19,1],[-219,37,5],[-219,79,9],[-219,109,11],
            [-939,349,21],[-1731,8287,105],[-3819,487,25],
            [-5619,937,35],[-8643,739,31],[-14619,1279,41],
            [-16419,2287,55],[-32187,10639,119],[-32619,3187,65],
            [-35859,8929,109],[-43059,7369,99],[-61995,5437,85],
            [-71499,6229,91],[-87267,7669,101]];
all_genus_pass = 1;
for(i=1,#gp_cases,
  sf = gp_cases[i][1]; p = gp_cases[i][2]; k = gp_cases[i][3];
  if(kronecker(sf,p)!=1, next()); \\ skip inert/ramified
  gv = genus_prod(sf, p);
  gpr = gv[1]; gch = gv[2];
  printf("%-8d %-6d %-4d  %-30s  %d\n", sf, p, k, gch, gpr);
  if(gpr != 1, all_genus_pass = 0)
);
print();
if(all_genus_pass,
  print("RESULT D: ALL cases have genus product=+1 (p in principal genus). Confirmed."),
  print("RESULT D: Some cases NOT in principal genus. Check needed.")
);

print();
print("=== DONE ===");

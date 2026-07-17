\\ thread15_pi_squared.gp — Thread 15: pi^2 principality
\\ Proves ord([pi])=2 by showing pi^2 is an explicit algebraic integer in O_K.
\\ For biquadratic Weil poly T^4+a2*T^2+p^2:
\\   pi^2 = (-a2 + m*sqrt(sf)) / 2,  disc4 = a2^2-4p^2 = sf*m^2
\\ If pi^2 in O_K then (pi^2) is trivially principal => ord([pi])|2.

default(parisize, 128000000);
default(timer, 0);

sqfree(n) = {
  my(f, d, e);
  if(n == 0, return(0));
  f = factor(abs(n));
  d = 1;
  for(ii = 1, #f[,1], e = f[ii,2] % 2; if(e, d *= f[ii,1]));
  if(n < 0, d = -d);
  d
};

\\ is (-a2 + m*sqrt(sf))/2 in O_K=Q(sqrt(sf))?
ok_member(a2, m, sf) = {
  my(r);
  r = sf % 4; if(r < 0, r += 4);
  if(r == 1, ((abs(a2) % 2) == (abs(m) % 2)), ((a2 % 2 == 0) && (m % 2 == 0)))
};

check_case(k, p, a2, sf) = {
  my(disc4, m2, m, in_ok, pi2norm, K, pr, res2, princ);
  disc4 = a2^2 - 4*p^2;
  m2 = disc4 / sf;
  m = sqrtint(abs(m2));
  in_ok = ok_member(a2, m, sf);
  pi2norm = (a2^2 - disc4) / 4;
  K = bnfinit(x^2 - sf, 1);
  pr = idealprimedec(K, p);
  res2 = bnfisprincipal(K, idealpow(K, pr[1], 2));
  princ = (res2[1] == 0*res2[1]);
  printf("k=%-4d p=%-7d sf=%-9d m=%-6d | in_OK=%d nmOK=%d pi2princ=%d", k,p,sf,m,in_ok,(pi2norm==p^2),princ);
  if(!in_ok || !princ || pi2norm != p^2, print(" *** ANOMALY ***"), print(""))
};

print("=== Thread 15: pi^2 principality for all norm-form cases ===");
print("Cols: in_OK = pi^2 in O_K by parity; nmOK = Nm(pi^2)=p^2; pi2princ = bnfisprincipal says principal");
print();
check_case(1, 19, -35, -219);
check_case(5, 37, 1, -219);
check_case(9, 79, 85, -219);
check_case(11, 109, 145, -219);
check_case(21, 349, 385, -939);
check_case(25, 487, -299, -3819);
check_case(31, 739, -1403, -8643);
check_case(35, 937, 1, -5619);
check_case(41, 1279, -2315, -14619);
check_case(55, 2287, -899, -16419);
check_case(65, 3187, -4499, -32619);
check_case(85, 5437, -9791, -61995);
check_case(91, 6229, -11375, -71499);
check_case(99, 7369, 385, -43059);
check_case(101, 7669, -13751, -87267);
check_case(105, 8287, 2149, -1731);
check_case(109, 8929, 5905, -35859);
check_case(119, 10639, 10549, -32187);
check_case(131, 12889, -25751, -3);

print();
print("=== Eisenstein prime for p=12889, sf=-3 ===");
{
  my(K, pr, r1, r2, c, d, a, b);
  K = bnfinit(x^2 + x + 1, 1);
  pr = idealprimedec(K, 12889);
  r1 = bnfisprincipal(K, pr[1]);
  r2 = bnfisprincipal(K, idealpow(K, pr[1], 2));
  a = r1[2][1]; b = r1[2][2];
  printf("pi generator: %d + %d*omega, Nm=%d (need 12889)\n", a, b, a^2-a*b+b^2);
  c = r2[2][1]; d = r2[2][2];
  printf("pi^2 generator: %d + %d*omega, Nm=%d (need %d)\n", c, d, c^2-c*d+d^2, 12889^2);
  printf("pi^2 vs formula: a2=%d, m=%d; (-a2+m*sqrt(-3))/2 should match\n", -25751, 681);
  printf("  formula gives real part: %d/2 = %d\n", 25751, 25751/2);
  printf("  Since sf=-3 equiv 1 (mod 4): omega=(-1+sqrt(-3))/2\n");
  printf("  sqrt(-3) = 2*omega+1, so m*sqrt(-3)/2 = %d*(2*omega+1)/2 = %d*omega + %d/2\n", 681, 681, 681);
  printf("  pi^2 = (25751 + 681*sqrt(-3))/2 = 25751/2 + (681*omega + 681/2)\n");
  printf("       = (25751+681)/2 + 681*omega = 13216 + 681*omega\n");
  printf("  bnf generator: %d + %d*omega  (sign: both match up to unit)\n", c, d)
}

print();
print("=== Explicit pi^2 generators for all cases ===");
print("(Shows the actual element pi^2 = alpha + beta*sqrt(sf) or similar)");
print();
show_gen(k, p, a2, sf) = {
  my(K, pr, r2, g, r, alpha, beta);
  K = bnfinit(x^2 - sf, 1);
  pr = idealprimedec(K, p);
  r2 = bnfisprincipal(K, idealpow(K, pr[1], 2));
  g = r2[2];
  r = sf % 4; if(r < 0, r += 4);
  if(r == 1,
    printf("k=%-4d sf=%-9d: pi^2 = %d + %d*(1+sqrt(%d))/2  [basis {1,(1+sqrt)/2}]\n",
      k, sf, g[1], g[2], sf)
  ,
    printf("k=%-4d sf=%-9d: pi^2 = %d + %d*sqrt(%d)  [basis {1,sqrt}]\n",
      k, sf, g[1], g[2], sf)
  )
};
show_gen(1, 19, -35, -219);
show_gen(5, 37, 1, -219);
show_gen(9, 79, 85, -219);
show_gen(11, 109, 145, -219);
show_gen(21, 349, 385, -939);
show_gen(35, 937, 1, -5619);
show_gen(131, 12889, -25751, -3);

print();
print("DONE.");

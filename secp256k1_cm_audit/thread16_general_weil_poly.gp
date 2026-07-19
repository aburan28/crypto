\\ thread16_general_weil_poly.gp
\\
\\ Thread 16: Verify order-2 Frobenius theorem extends to non-norm-form primes.
\\
\\ PRODUCT JACOBIAN CONSTRUCTION:
\\   E/F_p with Frobenius trace t => J = E x E^twist has Weil polynomial
\\     (T^2-t*T+p)(T^2+t*T+p) = T^4 + (2p-t^2)*T^2 + p^2
\\   so a2 = 2p - t^2.  D = t^2*(t^2-4p) < 0 by Hasse, so sf < 0.
\\   For t != 0: p does not divide a2, so Thread-15 proof (A)-(E) applies
\\   and predicts [P]^2 = 1 in Cl(Q(sqrt(sf))).
\\
\\ Run: gp -q thread16_general_weil_poly.gp

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

is_norm_form(p) = {
  my(val, q, sq);
  val = 4*p - 73;
  if(val <= 0 || val % 3 != 0, return(0));
  q = val \ 3;
  sq = sqrtint(q);
  (sq*sq == q) && (sq % 2 == 1)
};

\\ Verify [P]^2=1 for given prime p and biquadratic coefficient a2.
\\ Returns 1 if OK, 0 if not, -1 if degenerate.
check_order2(p, a2) = {
  my(D, sf, m2, m, K, Pp, P, P2, pr, ok);
  D  = a2^2 - 4*p^2;
  if(D == 0, return(-1));
  sf = sf_part(D);
  if(sf == 0, return(-1));
  m2 = D \ sf;
  if(m2 <= 0 || m2*sf != D, return(-1));
  m  = sqrtint(m2);
  if(m*m != m2, return(-1));
  if(a2 % p == 0, return(-1));

  K  = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  if(#Pp == 0, return(-1));

  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  pr = bnfisprincipal(K, P2);
  ok = (pr[1] == 0*pr[1]);
  [ok, sf, m, K.clgp.no, #Pp]
};

print("Thread 16: Order-2 Frobenius Ideal -- Extension to Non-Norm-Form Primes");
print("=========================================================================");
print("");
print("For each non-norm-form prime p, pick E/F_p with trace t != 0.");
print("Set a2 = 2p - t^2 (product Jacobian). Verify [P]^2=1 in Cl(Q(sqrt(sf))).");
print("");

{
  my(count, ok_count, p, found, a2_val, t_val, E, t, a2, res, ok, sf, m, h);
  count = 0;
  ok_count = 0;
  p = 50;

  while(count < 10,
    p = nextprime(p + 1);
    if(is_norm_form(p), next);

    \\ Find E/F_p with t != 0 using y^2 = x^3 + a*x + b
    found = 0;
    a2_val = 0;
    t_val  = 0;
    for(aa = 1, 50,
      if(found, break);
      for(bb = 1, 30,
        if(found, break);
        if((4*aa^3 + 27*bb^2) % p == 0, next);
        E = ellinit([aa, bb], p);
        t = ellap(E, p);
        if(t == 0, next);
        a2 = 2*p - t^2;
        if(a2 == 0 || a2 % p == 0, next);
        found = 1;
        a2_val = a2;
        t_val  = t
      )
    );

    if(!found, next);

    res = check_order2(p, a2_val);
    if(type(res) == "t_INT",
      printf("p=%-7d t=%-4d a2=%-8d: degenerate, skip\n", p, t_val, a2_val);
      next
    );

    ok = res[1];
    sf = res[2];
    m  = res[3];
    h  = res[4];

    printf("p=%-7d t=%-4d a2=%-8d sf=%-10d m=%-5d h=%-3d  [P]^2=1: %s\n",
      p, t_val, a2_val, sf, m, h, if(ok, "YES", "NO(!)"));

    count++;
    ok_count += ok
  );

  print("");
  printf("Tested %d non-norm-form product-Jacobian pairs.\n", count);
  printf("All [P]^2=1: %s\n", if(count == ok_count, "YES", "NO"));
  print("");
  print("NOTE: The proof (A)-(E) from Thread 15 is entirely general.");
  print("No norm-form condition enters the argument.");
}

print("");
print("DONE.");

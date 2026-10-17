crypto_class_certificate() = {
  my(d, form, n, fac, principal, candidate, prime);
  d = -6649310858224389603;
  form = Qfb(7,5,237475387793728201);
  n = 68162026;
  fac = factor(n);
  principal = Qfb(1,1,(1-d)/4);
  if(5^2-4*7*237475387793728201 != d, error("discriminant mismatch"));
  if(type(d)!="t_INT" || type(form)!="t_QFB", error("bad input type"));
  print("input_form=",form);
  print("discriminant=",d);
  print("order=",n);
  print("factor=",fac);
  print("principal=",principal);
  print("power_n=",qfbpow(form,n));
  if(qfbpow(form,n)!=principal, error("order does not annihilate class"));
  for(i=1,matsize(fac)[1],
    prime = fac[i,1];
    if(!isprime(prime),error("factor not prime"));
    candidate=qfbpow(form,n/prime);
    print("prime=",prime,"; prime_verified=1; power_n_over_prime=",candidate);
    if(candidate==principal,error("order not exact"));
  );
  print("exact_order_certificate=PASS");
};
iferr(crypto_class_certificate(), failure, print(Str(failure)); quit(1));
quit(0);

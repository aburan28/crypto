default(parisize, 64000000);
found_p = 0; found_b = 0; found_n = 0;
{
forprime(p = 5000, 200000,
  if(Mod(p,6) != 1, next);
  for(b = 1, 20,
    local(E, n);
    E = ellinit([0,b], p);
    n = ellcard(E);
    if(isprime(n) && n > 5000 && kronecker(-3, n) == 1,
      print("p=", p, " b=", b, " n=", n, " nbits=", #binary(n));
      if(found_p == 0, found_p = p; found_b = b; found_n = n)
    )
  )
);
}
print("Best: p=", found_p, " b=", found_b, " n=", found_n);

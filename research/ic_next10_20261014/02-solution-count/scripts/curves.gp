\\ Seeded search: for bits b in {20,24,28} and sign s in {+1,-1}, primes p = s mod 2^11 just above 2^b,
\\ three prime-order curves y^2 = x^3 + a x + b (one per prime, a,b from a deterministic sequence).
setrand(20261014);
out = [];
{
foreach([20,24,28], bits,
  foreach([1,-1], s,
    p = 2^bits + (s - 2^bits) % 2^11; if(p < 2^bits, p += 2^11);
    found = 0;
    while(found < 3,
      if(isprime(p),
        for(k=1, 400,
          a = random(p); b = random(p);
          if((4*a^3+27*b^2) % p == 0, next);
          E = ellinit([a,b], p); n = ellcard(E);
          if(isprime(n) && n != p,
            found++; out = concat(out, [[bits, s, p, a, b, n]]); break)));
      p += 2^11)));
}
foreach(out, c, printf("%d,%d,%d,%d,%d,%d\n", c[1], c[2], c[3], c[4], c[5], c[6]));

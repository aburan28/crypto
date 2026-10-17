\\ Validation curves: smallest prime p > 2^12 with p = s (mod 64), s = +1, -1, first prime-order curve.
setrand(20261014);
{
foreach([1,-1], s,
  p = 2^12 + (s - 2^12) % 64; while(!isprime(p), p += 64);
  for(k=1, 1000, a = random(p); b = random(p);
    if((4*a^3+27*b^2) % p == 0, next);
    n = ellcard(ellinit([a,b], p));
    if(isprime(n) && n != p, printf("12,%d,%d,%d,%d,%d\n", s, p, a, b, n); break)));
}

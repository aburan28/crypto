\\ Find (a,c) for correct Howe cover y^2 = x^6 + a*x^3 + c for j=0 curves.
\\ Tests small primes p == 1 mod 6 and b1 = 1..5.

default(parisize, 256000000);
default(timer, 0);

print("Finding correct Howe cover h=x^6+a*x^3+c for j=0 curves");
print("p   b1  t1   d   om     a     c      b2   b1+b2  b1*b2  a^2-4c");
print("--------------------------------------------------------------------");

results = List();

{
forprime(p = 7, 130,
  if(Mod(p,6) != 1, next);
  om_r = polrootsmod(x^2+x+1, p);
  if(#om_r == 0, next);
  om = lift(om_r[1]);
  d = 2;
  while(kronecker(d,p) != -1, d++);
  forstep(b1 = 1, 10, 1,
    E1 = ellinit([0, b1], p);
    n1 = ellcard(E1);
    t1 = p + 1 - n1;
    if(t1 == 0, next);
    coeff_x2 = 2*p - t1^2;
    expected_n = 1 + coeff_x2 + p^2;
    b2 = (d^3 * b1) % p;
    found_a = -1; found_c = -1;
    forstep(a = 0, p-1, 1,
      forstep(c = 1, p-1, 1,
        disc_c = poldisc(x^6 + a*x^3 + c) % p;
        if(disc_c == 0, next);
        h_c = Mod(x^6 + a*x^3 + c, p);
        cp = hyperellcharpoly(h_c);
        if(subst(cp, variable(cp), 1) == expected_n,
          found_a = a; found_c = c; break(2)
        )
      )
    );
    if(found_a >= 0,
      disc_ac = (found_a^2 - 4*found_c) % p;
      printf("%3d  %2d  %3d  %2d  %3d   %4d  %4d   %4d   %4d  %5d  %5d\n",
             p, b1, t1, d, om, found_a, found_c, b2, (b1+b2)%p, (b1*b2)%p, disc_ac);
      listput(results, [p, b1, t1, d, om, found_a, found_c, b2])
    )
  )
)
}

print("");
print("Pattern check: is c = b1*b2 mod p?   is a = b1+b2 mod p?");
print("or: is c = sqrt(b1*b2)?  is disc(a^2-4c) related to -3?");
print("");
print("Per-case analysis:");
for(i = 1, #results,
  r = results[i];
  p_r=r[1]; b1_r=r[2]; t1_r=r[3]; d_r=r[4]; om_r=r[5]; a_r=r[6]; c_r=r[7]; b2_r=r[8];
  \\ Check: does disc(a^2-4c) relate to b1*d^3?
  disc = Mod(a_r^2 - 4*c_r, p_r);
  \\ is disc a square?
  is_sq = (kronecker(lift(disc), p_r) == 1);
  \\ ratio c_r / b2_r
  c_over_b2 = lift(Mod(c_r, p_r) / Mod(b2_r, p_r));
  \\ ratio a / (b1+b2)
  b1_plus_b2 = (b1_r + b2_r) % p_r;
  a_over_sum = if(b1_plus_b2 != 0, lift(Mod(a_r, p_r)/Mod(b1_plus_b2, p_r)), -1);
  printf("p=%3d b1=%2d b2=%3d a=%3d c=%3d  c/b2=%3d  a/(b1+b2)=%d  disc_sq=%d\n",
         p_r, b1_r, b2_r, a_r, c_r, c_over_b2, a_over_sum, is_sq)
);

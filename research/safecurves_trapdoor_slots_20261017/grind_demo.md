# Grind demonstration: strong quadratic twist, weak sextic sibling

Executed by `sextic_grind_demo.py`. For each size, primes p = 1 mod 6 were
scanned (all six sextic classes of y^2 = x^3 + b over F_p computed by SEA);
the first (p, class) pair is reported where the published class and its
quadratic twist both have largest prime factor >= 0.79 * bitlen(p) bits while
some other sextic class has largest prime factor <= 0.40 * bitlen(p).

## 48 bits (found after 2 primes = 1 mod 6)

- p = 281474976710731
- publish b in class d = 3: order 281474959933501, largest factor 48 bits
- its quadratic twist (class pair): order 281474993487963, largest factor 40 bits
- weak sextic sibling d = 2: order 281474959933531, largest factor 13 bits -> rho cost ~ 2^6.5

| class d | order | largest factor (bits) |
|---:|---|---:|
| 1 | 281474943156300 | 20 |
| 3 | 281474959933501 | 48 |
| 2 | 281474959933531 | 13 |
| 9 | 281474993487933 | 24 |
| 4 | 281474993487963 | 40 |
| 8 | 281475010265164 | 47 |

## 96 bits (found after 3 primes = 1 mod 6)

- p = 79228162514264337593543950963
- publish b in class d = 4: order 79228162514263994975295378741, largest factor 95 bits
- its quadratic twist (class pair): order 79228162514264680211792523187, largest factor 79 bits
- weak sextic sibling d = 7: order 79228162514264122064149793307, largest factor 28 bits -> rho cost ~ 2^14.0

| class d | order | largest factor (bits) |
|---:|---|---:|
| 8 | 79228162514263779445901221084 | 42 |
| 4 | 79228162514263994975295378741 | 95 |
| 7 | 79228162514264122064149793307 | 28 |
| 2 | 79228162514264553122938108621 | 68 |
| 5 | 79228162514264680211792523187 | 79 |
| 1 | 79228162514264895741186680844 | 31 |


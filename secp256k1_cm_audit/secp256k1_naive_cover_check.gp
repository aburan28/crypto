\\ Quick check: is y^2=(x^3+7)(x^3+189) the correct Howe-glued cover for secp256k1?
\\ secp256k1: y^2 = x^3 + 7 over F_p where p = 2^256 - 2^32 - 977.

default(parisize, 256000000); default(timer, 0);

print("================================================================");
print("secp256k1 naive cover check: y^2=(x^3+7)(x^3+189)");
print("================================================================");
print("");

p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
b = 7; b2 = 189;

E1 = ellinit([0, b], p_secp);
t1 = p_secp + 1 - ellcard(E1);
print("E1: y^2 = x^3 + 7 over F_p_secp");
print("trace of Frobenius t1 = ", t1);
print("");

E2 = ellinit([0, b2], p_secp);
t2 = p_secp + 1 - ellcard(E2);
print("E2: y^2 = x^3 + 189 over F_p_secp");
print("trace of Frobenius t2 = ", t2);
print("Is E2 the quadratic twist of E1? (t2 == -t1): ", t2 == -t1);
print("");

\\ Is -b = -7 a cube mod p_secp?
neg_b = p_secp - b;
cube_test = lift(Mod(neg_b, p_secp)^((p_secp-1)/3));
print("Is -b = p-7 a cube mod p_secp? (p-1)/3 power = 1 means yes.");
print("(-b)^((p-1)/3) mod p = ", cube_test, "  (1 = cube, else not)");
print("");

\\ Can we compute char poly of the genus-2 curve?
\\ This is too large for hyperellcharpoly to verify directly (p is 256-bit).
\\ But we can check the Weil poly formula: for Jac ~ E1 x E2,
\\ the Weil poly = (x^2 - t1*x + p)(x^2 - t2*x + p) = (x^2-t1*x+p)(x^2+t1*x+p)
\\ = x^4 + (2*p - t1^2)*x^2 + p^2
weil_mid = 2*p_secp - t1^2;
print("Expected Weil poly if Jac ~ E1 x E2:");
print("  x^4 + (2p - t^2)*x^2 + p^2");
print("  2p - t^2 = ", weil_mid);
print("  (sign of 2p-t^2: ", if(weil_mid > 0, "positive", "negative"), ")");
print("");

\\ For the naive cover h = (x^3+7)(x^3+189), verify it's smooth (disc != 0):
h_secp = (x^3 + b) * (x^3 + b2);
J10 = poldisc(h_secp);
print("Discriminant of h = (x^3+7)(x^3+189): disc = ", J10);
print("  => smooth genus-2 curve over Z: ", J10 != 0);
print("");

\\ The key missing check: does hyperellcharpoly(h_secp mod p_secp)
\\ equal x^4 + weil_mid*x^2 + p_secp^2?
\\ This is expensive for 256-bit p. Instead, verify a SMALLER prime
\\ p' ≡ p_secp (mod relevant structure) and b'=7, b2'=189.

\\ Find a small prime p' ≡ 1 (mod 6) where E: y^2=x^3+7 has the same
\\ qualitative properties (non-rational 2-torsion, quadratic twist is y^2=x^3+189).
print("---- Proxy test over small prime p' ≡ 1 mod 6, b=7, b2=189 ----");
{
  forprime(pp = 100, 2000,
    if(pp % 6 != 1, next);
    E_small = ellinit([0, 7], pp);
    t_small = pp + 1 - ellcard(E_small);
    E2_small = ellinit([0, 189], pp);
    t2_small = pp + 1 - ellcard(E2_small);
    if(t2_small != -t_small, next);   \\ require quadratic twist
    if(t_small == 0, next);           \\ skip supersingular
    \\ Check the naive cover char poly
    h_small = Mod(1,pp)*x^6 + Mod(196,pp)*x^3 + Mod(1323,pp);
    if(poldisc(h_small) == 0, next);
    cp = hyperellcharpoly(h_small);
    expected_mid = 2*pp - t_small^2;
    actual_mid = polcoeff(cp, 2);
    match = (lift(actual_mid) == expected_mid % pp);
    print("  p'=", pp, "  t=", t_small, "  expected_mid=", expected_mid, "  actual_mid=", lift(actual_mid), "  match=", match);
    if(match, break)
  )
}

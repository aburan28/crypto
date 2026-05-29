\\ Flat (no nested loops) check of Rosenhain applicability
\\ Hardcode a (p,b) where (-b) IS a cubic residue to verify formula.
\\
\\ Pre-computed: p=13, b=5: (-5 mod 13 = 8), 8^{(13-1)/3} = 8^4 = 4096 mod 13 = 4
\\   Not 1. Try p=13,b=1: (-1 mod 13=12), 12^4=20736 mod 13. 12=-1, (-1)^4=1. OK!
\\
\\ p=13, b=1: y^2 = x^3 + 1 (j=0, p=13=1 mod 3).
\\ x^3 = -1 mod 13: x^3 = 12 mod 13.
\\ 12 = -1, so x = -1 = 12 (since 12^3 = -1 = 12 mod 13), and...
\\ Cube roots of -1 mod 13: 12, 12*om, 12*om^2 where om^3=1.
\\ om satisfies om^2+om+1=0 mod 13. Roots of x^2+x+1 mod 13:
\\ x = (-1 pm sqrt(-3))/2 mod 13. -3 mod 13 = 10. sqrt(10) mod 13?
\\ 6^2=36=10 mod 13. So om = (-1+6)/2 = 5/2. 2^{-1} mod 13 = 7. om = 5*7=35=9 mod 13.
\\ Check: 9^2+9+1 = 81+9+1 = 91 = 7*13 = 0 mod 13. Correct.
\\ Cube roots of -1: 12, 12*9=108=4 mod 13, 12*9^2=12*81=12*3=36=10 mod 13.
\\ Check: 4^3=64=12 mod 13. 10^3=1000=1000-76*13=1000-988=12 mod 13. Both are -1=12. Correct.
\\
\\ So for (p=13, b=1): 2-torsion x-coords = {12, 4, 10}. All F_13-rational.
\\ omega = 9 (cube root of unity mod 13).
\\
\\ Now take d=2 (non-square mod 13? 2^6 mod 13 = 64 mod 13 = 12 = -1. Yes, 2 is a non-square).
\\ E1: y^2 = x^3 + 1 mod 13
\\ E2: y^2 = x^3 + d^3*1 = x^3 + 8 mod 13
\\
\\ Run: gp -q chlrs_flat_check.gp

default(parisize, 32000000);

p = 13;
b = 1;
d_ns = 2;
om = 9;

print("=== Rosenhain validity check for (p=13, b=1, d=2, om=9) ===");
printf("p=%d, b=%d (j=0), d=%d (non-square)\n", p, b, d_ns);
printf("omega = %d, omega^2+omega+1 = %d mod p\n", om, lift(Mod(om^2+om+1,p)));
printf("(-b)^{(p-1)/3} = %d (should be 1 for cubic res)\n", lift(Mod(-b,p)^((p-1)/3)));

rts = polrootsmod(x^3 + b, p);
printf("Roots of x^3+1 mod 13: %d roots\n", #rts);
printf("  %d, %d, %d\n", lift(rts[1]), lift(rts[2]), lift(rts[3]));

E1 = ellinit([0, b], p);
b2 = lift(Mod(d_ns^3, p));
E2 = ellinit([0, b2], p);
t1 = p + 1 - ellcard(E1);
t2 = p + 1 - ellcard(E2);
printf("b2 = d^3*b mod p = %d\n", b2);
printf("t(E1) = %d, t(E2) = %d (expected %d)\n", t1, t2, -t1);

d = Mod(d_ns, p);
w = Mod(om, p);
l1 = (w^2-1)*(d-w) / (w*(w-1)*(d-1));
l2 = (d*w-1)*(d-w) / (w*(d-1)^2);
l3 = (d*w^2-1)*(d-w) / (w*(d*w-1)*(d-1));
printf("Rosenhain l1=%d, l2=%d, l3=%d\n", lift(l1), lift(l2), lift(l3));

h = Mod(1,p)*x*(x-1)*(x-l1)*(x-l2)*(x-l3);
cp = hyperellcharpoly(h);
n_jac = subst(cp, variable(cp), 1);
target_n = (p+1-t1)*(p+1+t1);
printf("Frobenius char poly: %Ps\n", cp);
printf("#Jac(C) = %d, target = %d, match = %d\n", n_jac, target_n, n_jac == target_n);
print("");

\\ secp256k1 final verdict
print("=== secp256k1 verdict ===");
p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
printf("(-7)^{(p-1)/3} mod p_secp = 1? %d\n", lift(Mod(-7,p_secp)^((p_secp-1)/3)) == 1);
printf("x^3+7 roots mod p_secp: %d\n", #polrootsmod(x^3+7, p_secp));
print("CONCLUSION: 2-torsion of secp256k1 NOT F_p-rational.");
print("  F_p Rosenhain formula CANNOT be applied directly.");
print("  Required: either F_{p^3} Rosenhain or Mestre reconstruction.");

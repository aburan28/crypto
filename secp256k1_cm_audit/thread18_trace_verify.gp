\\ Thread 18: Verify correct CM trace assignment for secp256k1 sextic twists.
\\ Run: gp -q thread18_trace_verify.gp

default(parisize, 1024000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t_secp = p + 1 - n_secp;
s = sqrtint((4*p - t_secp^2)/3);

g_gen = znprimroot(p);
u = lift(g_gen ^ ((p - 1) / 6));

\\ Six b-values
b0 = 7;
b1 = lift(Mod(7,p) * Mod(u,p));
b2 = lift(Mod(7,p) * Mod(u,p)^2);
b3 = p - 7;
b4 = lift(Mod(7,p) * Mod(u,p)^4);
b5 = lift(Mod(7,p) * Mod(u,p)^5);

\\ Six abstract CM traces
T1 = t_secp;
T2 = (t_secp - 3*s) / 2;
T3 = -(t_secp + 3*s) / 2;
T4 = -t_secp;
T5 = (3*s - t_secp) / 2;
T6 = (t_secp + 3*s) / 2;

print("=== CM trace candidates ===");
print("T1 = t   = ", T1, "  (N1 mod 4 = ", (p+1-T1)%4, ")");
print("T2 =     = ", T2, "  (N2 mod 4 = ", (p+1-T2)%4, ")");
print("T3 =     = ", T3, "  (N3 mod 4 = ", (p+1-T3)%4, ")");
print("T4 = -t  = ", T4, "  (N4 mod 4 = ", (p+1-T4)%4, ")");
print("T5 =     = ", T5, "  (N5 mod 4 = ", (p+1-T5)%4, ")");
print("T6 =     = ", T6, "  (N6 mod 4 = ", (p+1-T6)%4, ")");
print("");
print("k=1 and k=4 have [1,1,1] 2-torsion so 4|N.");
print("N with 4|N: T3 (N mod 4=", (p+1-T3)%4, ") and T6 (N mod 4=", (p+1-T6)%4, ")");
print("");

\\ Find a point on E_{k=1}: y^2 = x^3 + b1
\\ Try small x values
Ek1 = ellinit([Mod(0,p), Mod(b1,p)]);
xtest = 2;
rhs = lift(Mod(xtest,p)^3 + Mod(b1,p));
while (kronecker(rhs, p) != 1,
    xtest = xtest + 1;
    rhs = lift(Mod(xtest,p)^3 + Mod(b1,p))
);
ytest = lift(sqrtmod(rhs, p));
Pk1 = [Mod(xtest,p), Mod(ytest,p)];
print("Point on E_{k=1}: x=", xtest, " y=", ytest);

N_T3 = p + 1 - T3;
N_T6 = p + 1 - T6;
res3 = ellmul(Ek1, Pk1, N_T3);
res6 = ellmul(Ek1, Pk1, N_T6);
print("N(T3)*P on E_{k=1}: x-coord mod 10^10 = ", res3[1] % 10^10);
print("N(T6)*P on E_{k=1}: x-coord mod 10^10 = ", res6[1] % 10^10);
print("Is N(T3)*P = O? ", (res3 == [0]) || (type(res3[1]) == "t_INT" && res3[1] == 0));
print("Is N(T6)*P = O? ", (res6 == [0]) || (type(res6[1]) == "t_INT" && res6[1] == 0));

\\ More reliable: check if result is zero (point at infinity in affine coords is [0])
print("");
print("Test: N(T3)*P:");
print(res3);
print("Test: N(T6)*P:");
print(res6);
print("");

\\ Find a point on E_{k=2}: y^2 = x^3 + b2
Ek2 = ellinit([Mod(0,p), Mod(b2,p)]);
xtest2 = 2;
rhs2 = lift(Mod(xtest2,p)^3 + Mod(b2,p));
while (kronecker(rhs2, p) != 1,
    xtest2 = xtest2 + 1;
    rhs2 = lift(Mod(xtest2,p)^3 + Mod(b2,p))
);
ytest2 = lift(sqrtmod(rhs2, p));
Pk2 = [Mod(xtest2,p), Mod(ytest2,p)];
print("Point on E_{k=2}: x=", xtest2, " y=", ytest2);

N_T2 = p + 1 - T2;
N_T5 = p + 1 - T5;
res_T2 = ellmul(Ek2, Pk2, N_T2);
res_T5 = ellmul(Ek2, Pk2, N_T5);
print("Test: N(T2)*P on E_{k=2}:");
print(res_T2);
print("Test: N(T5)*P on E_{k=2}:");
print(res_T5);
print("");

print("=== CORRECT order assignments ===");
print("(From scalar-mult tests above: the O result identifies the correct trace)");

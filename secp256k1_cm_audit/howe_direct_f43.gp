\\ howe_direct_f43.gp  (PARI 2.15.4 compat: all logic in functions, no multiline for-bodies)
\\ Run: gp -q howe_direct_f43.gp

default(parisize, 512000000); default(timer, 0);
p = 43;

\\ ------------------------------------------------------------------
\\ Helper: Legendre symbol (a|p) via Euler criterion
\\ ------------------------------------------------------------------
leg(a) = lift(Mod(a,p)^((p-1)/2));  \\ 1=square, 42=-1=non-square

\\ ------------------------------------------------------------------
\\ Helper: count affine pts on C: y^2=(x^3+7)(x^3+13) over F_p
\\ ------------------------------------------------------------------
count_affine_Fp(px) = {
    my(v, r);
    v = lift(Mod((px^3+7)*(px^3+13), p));
    if (v==0, r=1, if(leg(v)==1, r=2, r=0));
    r
};

count_C_Fp() = {
    my(total, xx);
    total = 0;
    for (xx=0, p-1, total = total + count_affine_Fp(xx));
    total + 2   \\ +2 pts at infinity (monic degree 6)
};

\\ ------------------------------------------------------------------
\\ Helper: count affine pts over F_{43^2} = F_43[i], i^2=2
\\ Element a+bi; x^3 = (a^3+6ab^2) + (3a^2b+2b^3)*i
\\ is (pre+pim*i) a square in F_{43^2}?
\\   -> z^{(p^2-1)/2} = 1 in F_{43^2}
\\ ------------------------------------------------------------------
is_square_Fp2(pre, pim) = {
    my(z, zp);
    if (pre==0 && pim==0, return(1));     \\ zero is a square
    z = Mod(pre + pim*x, Mod(x^2-2, p));
    zp = lift(lift(z^((p^2-1)/2)));
    zp == 1
};

count_affine_Fp2(aa, bb) = {
    my(re3, im3, re7, im7, re13, im13, pre, pim);
    re3  = lift(Mod(aa^3 + 6*aa*bb^2, p));
    im3  = lift(Mod(3*aa^2*bb + 2*bb^3, p));
    re7  = lift(Mod(re3+7, p));  im7  = im3;
    re13 = lift(Mod(re3+13,p));  im13 = im3;
    pre  = lift(Mod(re7*re13 + 2*im7*im13, p));
    pim  = lift(Mod(re7*im13 + re13*im7, p));
    my(v);
    if (pre==0 && pim==0, v=1, if(is_square_Fp2(pre,pim), v=2, v=0));
    v
};

count_C_Fp2() = {
    my(total, aa, bb);
    total = 0;
    for (aa=0, p-1, for (bb=0, p-1, total = total + count_affine_Fp2(aa, bb)));
    total + 2
};

\\ ==================================================================
print("=======================================================");
print("Howe direct: y^2=(x^3+7)(x^3+13) = x^6+20x^3+5 / F_43");
print("=======================================================");
print("");

\\ --- Sub-factor identification ---
a_h=20; b_h=5;
disc_sf = lift(Mod(a_h^2 - 4*b_h, p));
sq_sf   = lift(Mod(disc_sf,p)^((p+1)/4));
inv2    = lift(Mod(2,p)^(-1));
r1 = lift(Mod((-a_h + sq_sf)*inv2, p));
r2 = lift(Mod((-a_h - sq_sf)*inv2, p));
mr1 = lift(Mod(-r1,p)); mr2 = lift(Mod(-r2,p));
t_r1 = p+1 - ellcard(ellinit([0,mr1],p));
t_r2 = p+1 - ellcard(ellinit([0,mr2],p));
t1   = p+1 - ellcard(ellinit([0,7],p));
t2   = p+1 - ellcard(ellinit([0,13],p));

print("Sub-factors of y^2=x^6+20x^3+5:");
print("  E_{r1}: y^2=x^3+", mr1, "  trace=", t_r1);
print("  E_{r2}: y^2=x^3+", mr2, "  trace=", t_r2);
print("Target E1: y^2=x^3+7  trace=", t1, "  |E1|=", p+1-t1);
print("Target E2: y^2=x^3+13 trace=", t2, "  |E2|=", p+1-t2);
sf_match = (Set([mr1,mr2])==Set([7,13]));
print("Sub-factors match {E1,E2}? ", sf_match);
print("");

\\ --- Point count C(F_43) ---
print("Counting C(F_43)...");
cnt1 = count_C_Fp();
b1   = p+1 - cnt1;
print("#C(F_43) = ", cnt1, "  b1=", b1, "  (expected 44, b1=0)");
print("");

\\ --- Point count C(F_{43^2}) ---
print("Counting C(F_{43^2}) [~1849 iterations, may take ~30s]...");
cnt2 = count_C_Fp2();
print("#C(F_{43^2}) = ", cnt2);

\\ Recover b2 from Newton's identity: sum_i alpha_i^2 = b1^2 - 2*b2
\\ #C(F_p^2) = p^2 + 1 - (b1^2 - 2*b2)  =>  b2 = (b1^2 - (p^2+1-cnt2)) / 2
p2_psum = p^2 + 1 - cnt2;
b2 = (b1^2 - p2_psum) / 2;
print("b2 = ", b2, "  (expected -83)");
print("");

\\ Jacobian point count
Jac_Fp = 1 - b1 + b2 - b1*p + p^2;   \\ L(1) = |Jac(C)(F_p)|
print("|Jac(C)(F_43)| = L(1) = ", Jac_Fp);
print("Expected |E1|*|E2| = ", (p+1-t1)*(p+1-t2));
print("L-poly = L_E1*L_E2 (b1=0, b2=-83)? ", b1==0 && b2==-83);
print("");

\\ --- Igusa J2, J10 ---
\\ f = x^6+20x^3+5: a6=1,a5=0,a4=0,a3=20,a2=0,a1=0,a0=5
J2  = lift(Mod(-120*5*1 + 3*20^2, p));
J10 = lift(Mod(poldisc(Pol([1,0,0,20,0,0,5])), p));
print("Igusa J2  mod 43 = ", J2);
print("Igusa J10 mod 43 = ", J10, "  (nonzero => smooth)");
print("");

print("=======================================================");
print("SUMMARY:");
print("  Sub-factors match (E1,E2)? ", sf_match);
print("  #C(F_43)   = ", cnt1, "  b1=", b1);
print("  #C(F_43^2) = ", cnt2, "  b2=", b2);
print("  |Jac(C)(F_43)| = ", Jac_Fp, " vs |E1||E2|=", (p+1-t1)*(p+1-t2));
print("  L-poly = L_E1*L_E2? ", b1==0 && b2==-83);
print("  J2=", J2, "  J10=", J10);
print("=======================================================");

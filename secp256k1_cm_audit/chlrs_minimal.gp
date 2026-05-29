\\ Minimal Rosenhain formula verification
\\ p=13, b=1, d=2, omega=9: all conditions hold (cubic residue OK)
\\ Run: gp -q chlrs_minimal.gp

default(parisize, 32000000);

p = 13;

\\ 2-torsion x-coords of y^2 = x^3 + 1 mod 13: roots of x^3 = -1 = 12 mod 13
\\ Roots: 12, 4, 10 (verified: 12^3=1728=12*144=12*11*13+..., 12^3 mod 13 = 12. Wait:
\\ 12 = -1, (-1)^3 = -1 = 12 mod 13. Yes. 4^3 = 64 = 64-4*13 = 12. 10^3=1000, 1000 mod 13:
\\ 13*76=988, 1000-988=12. All correct.)

print("=== Verify Rosenhain formula for p=13, b=1, d=2 ===");
print("Curve: y^2 = x^3 + 1 mod 13, trace = ?");
E1 = ellinit([0, 1], 13);
t1 = 13 + 1 - ellcard(E1);
print("trace t1 = ", t1);

E2 = ellinit([0, 8], 13);
t2 = 13 + 1 - ellcard(E2);
print("trace t2 = ", t2, " (expected ", -t1, ")");
print("target #Jac = (13+1-t1)*(13+1+t1) = ", (14-t1)*(14+t1));

\\ Rosenhain invariants with omega=9, d=2
\\ l1 = 11, l2 = 7, l3 = 3 (computed manually and by script)
l1 = 11; l2 = 7; l3 = 3;

\\ Construct h as polynomial over Z, then pass p to hyperellcharpoly
h_int = x*(x-1)*(x-l1)*(x-l2)*(x-l3);
print("h(x) over Z = ", h_int);

\\ Method 1: pass polynomial over Z with p as second arg
cp1 = hyperellcharpoly(h_int, p);
nj1 = subst(cp1, variable(cp1), 1);
print("Method 1 (Z poly + explicit p): char poly = ", cp1);
print("  #Jac = ", nj1, "  target = ", (14-t1)*(14+t1), "  match = ", nj1 == (14-t1)*(14+t1));

\\ Method 2: reduce coefficients mod p, pass as Fp polynomial
h_fp = Mod(1, p)*x*(x-1)*(x-l1)*(x-l2)*(x-l3);
cp2 = hyperellcharpoly(h_fp);
nj2 = subst(cp2, variable(cp2), 1);
print("Method 2 (Fp poly, no explicit p): char poly = ", cp2);
print("  #Jac = ", nj2, "  target = ", (14-t1)*(14+t1), "  match = ", nj2 == (14-t1)*(14+t1));

\\ Also try the y^2 = DEGREE 6 form with all 6 branch points
\\ Branch points from E1: 12, 4, 10 (x^3 = -1 = 12 mod 13)
\\ Branch points from E2: 2*alpha where alpha^3 = -1, i.e., 2*12=24=11, 2*4=8, 2*10=20=7 mod 13
\\ Wait: E2: y^2 = x^3 + 8, 2-torsion: x^3 = -8 = 5 mod 13.
\\ -8 = -8 + 13 = 5. Roots of x^3 = 5 mod 13:
rts2 = polrootsmod(x^3 - 5, p);
print("2-torsion of E2 (roots of x^3 = 5 mod 13): ", rts2);

\\ So E2 branch points: distinct from E1 branch points {12,4,10}.
\\ The full sextic (degree 6) model:
\\ y^2 = (x-12)(x-4)(x-10) * (product of E2 branch points)
if (#rts2 == 3,
  r2_1 = lift(rts2[1]);
  r2_2 = lift(rts2[2]);
  r2_3 = lift(rts2[3]);
  print("E2 2-torsion: ", r2_1, ", ", r2_2, ", ", r2_3);
  h6 = (x-12)*(x-4)*(x-10)*(x-r2_1)*(x-r2_2)*(x-r2_3);
  print("Degree-6 model: y^2 = ", h6, " mod ", p);
  cp6 = hyperellcharpoly(h6, p);
  nj6 = subst(cp6, variable(cp6), 1);
  print("Char poly (degree-6 model): ", cp6);
  print("  #Jac = ", nj6, "  target = ", (14-t1)*(14+t1), "  match = ", nj6 == (14-t1)*(14+t1));
,
  print("E2 2-torsion not fully F_p-rational (#rts2 = ", #rts2, ")");
);

print("");
print("=== secp256k1 2-torsion verdict ===");
p_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
rts_secp = polrootsmod(x^3 + 7, p_secp);
print("Roots of x^3+7 mod p_secp: ", #rts_secp, " (0 = 2-torsion not F_p-rational)");
print("BLOCKED: secp256k1 Rosenhain formula requires F_{p^3} or Mestre approach.");

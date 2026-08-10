\\ ============================================================
\\ Mestre's Step 2: reconstruct a genus-2 curve from Igusa-Clebsch
\\ invariants (I2, I4, I6, I10) over Q.
\\ ============================================================
\\
\\ Closes gap #3 from RESEARCH_MESTRE_HOWE.md §7 ("Medium value,
\\ high cost: Mestre's Step 2 ... ~500 lines of PARI; week-long
\\ effort"): given (I2,I4,I6,I10), build Mestre's conic, find a
\\ rational point (PARI's qfsolve), parametrize it (PARI's
\\ qfparam), and assemble the reconstructed sextic f(x) with
\\ y^2 = f(x) having the given invariants (up to the standard
\\ weighted scaling (lambda^2,lambda^4,lambda^6,lambda^10)).
\\
\\ Source of the formulas: SageMath
\\ sage/schemes/hyperelliptic_curves/mestre.py
\\ (Mestre_conic + HyperellipticCurve_from_invariants), which in
\\ turn cites Lauter-Yang 2001 p.956-957 and Mestre 1991.  Fetched
\\ verbatim from
\\ https://raw.githubusercontent.com/sagemath/sage/develop/src/sage/schemes/hyperelliptic_curves/mestre.py
\\ on 2026-08-10 (arxiv.org / eprint.iacr.org are blocked by this
\\ environment's egress proxy; raw.githubusercontent.com is not).
\\
\\ This script implements the algorithm ONLY over Q (and number
\\ fields, via qfsolve's Hasse-Minkowski machinery). It does NOT
\\ implement the finite-field case: PARI's qfsolve requires an
\\ integer matrix (confirmed by direct test: qfsolve on a Mod(.,p)
\\ matrix raises "incorrect type in qfsolve [integer matrix]").
\\ Point-finding on a ternary conic over F_p needs a different,
\\ from-scratch algorithm (diagonalize the form mod p by completing
\\ the square, then find an isotropic vector — always possible for
\\ p odd since every ternary quadratic form over a finite field is
\\ isotropic). That is the concrete next step; see the bottom of
\\ this file and RESEARCH_MESTRE_HOWE.md.
\\
\\ Run: gp -q mestre_reconstruction.gp > mestre_reconstruction_output.txt
\\ (must run from this directory: uses read("igusa_clebsch.gp"))

default(parisize, 256000000);
default(timer, 0);

read("igusa_clebsch.gp");
read("igusa_clebsch_complete.gp");

\\ ------------------------------------------------------------
\\ Step 1 of Mestre_conic: (I2,I4,I6,I10) -> (x,y,z)
\\ ------------------------------------------------------------
Mestre_xyz(I2,I4,I6,I10) = {
    my(xx,yy,zz);
    xx = 8*(1 + 20*I4/I2^2)/225;
    yy = 16*(1 + 80*I4/I2^2 - 600*I6/I2^3)/3375;
    zz = -64*(-10800000*I10/I2^5 - 9 - 700*I4/I2^2 + 3600*I6/I2^3
             + 12400*I4^2/I2^4 - 48000*I4*I6/I2^5)/253125;
    [xx,yy,zz]
};

\\ Mestre's conic as a symmetric 3x3 Gram matrix L: the conic is
\\ Q(u,v,w) = u~ L u (PARI's qfsolve/qfparam convention already
\\ double-counts off-diagonal entries, matching Sage's Conic(k,L)).
Mestre_L(xx,yy,zz) = {
    matrix(3,3,i,j,
      if(i==1&&j==1, xx+6*yy,
      if((i==1&&j==2)||(i==2&&j==1), 6*xx^2+2*yy,
      if((i==1&&j==3)||(i==3&&j==1), 2*zz,
      if(i==2&&j==2, 2*zz,
      if((i==2&&j==3)||(i==3&&j==2), 9*xx^3+4*xx*yy+6*yy^2,
      if(i==3&&j==3, 6*xx^2*yy+2*yy^2+3*xx*zz, 0)))))))
};

\\ ------------------------------------------------------------
\\ Step 4 of Mestre's algorithm: assemble the sextic from the
\\ conic parametrization F1,F2,F3 (each a quadratic in one
\\ variable) via the c_ijk rational functions of (x,y,z).
\\ ------------------------------------------------------------
Mestre_sextic_from_F(xx,yy,zz,F1,F2,F3) = {
    my(c111,c112,c113,c122,c123,c133,c222,c223,c233,c333);
    c111 = 12*xx*yy - 2*yy/3 - 4*zz;
    c112 = -18*xx^3 - 12*xx*yy - 36*yy^2 - 2*zz;
    c113 = -9*xx^3 - 36*xx^2*yy - 4*xx*yy - 6*xx*zz - 18*yy^2;
    c122 = c113;
    c123 = -54*xx^4 - 36*xx^2*yy - 36*xx*yy^2 - 6*xx*zz - 4*yy^2 - 24*yy*zz;
    c133 = -27*xx^4/2 - 72*xx^3*yy - 6*xx^2*yy - 9*xx^2*zz - 39*xx*yy^2 - 36*yy^3 - 2*yy*zz;
    c222 = -27*xx^4 - 18*xx^2*yy - 6*xx*yy^2 - 8*yy^2/3 + 2*yy*zz;
    c223 = 9*xx^3*yy - 27*xx^2*zz + 6*xx*yy^2 + 18*yy^3 - 8*yy*zz;
    c233 = -81*xx^5/2 - 27*xx^3*yy - 9*xx^2*yy^2 - 4*xx*yy^2 + 3*xx*yy*zz - 6*zz^2;
    c333 = 27*xx^4*yy/2 - 27*xx^3*zz/2 + 9*xx^2*yy^2 + 3*xx*yy^3 - 6*xx*yy*zz + 4*yy^3/3 - 10*yy^2*zz;
    c111*F1^3 + c112*F1^2*F2 + c113*F1^2*F3 + c122*F1*F2^2 + c123*F1*F2*F3
      + c133*F1*F3^2 + c222*F2^3 + c223*F2^2*F3 + c233*F2*F3^2 + c333*F3^3
};

\\ Full pipeline over Q. Returns a sextic in variable 'tvar', or
\\ raises an error (with the local obstruction prime, per PARI's
\\ qfsolve convention) if no rational point exists on the conic.
Mestre_reconstruct(I2,I4,I6,I10) = {
    my(xyz,xx,yy,zz,L,sol,par,F1,F2,F3,f,t);
    xyz = Mestre_xyz(I2,I4,I6,I10);
    xx=xyz[1]; yy=xyz[2]; zz=xyz[3];
    L = Mestre_L(xx,yy,zz);
    sol = qfsolve(L);
    if(type(sol) == "t_INT",
        error("No rational point on Mestre conic (local obstruction at ", sol, ")"));
    par = qfparam(L, sol);
    t = 'tvar;
    F1 = par[1,1]*t^2 + par[1,2]*t + par[1,3];
    F2 = par[2,1]*t^2 + par[2,2]*t + par[2,3];
    F3 = par[3,1]*t^2 + par[3,2]*t + par[3,3];
    f = Mestre_sextic_from_F(xx,yy,zz,F1,F2,F3);
    f * denominator(f)
};

\\ ============================================================
\\ Tests
\\ ============================================================

print("================================================================");
print("Mestre reconstruction: (I2,I4,I6,I10) -> genus-2 curve over Q");
print("================================================================");
print("");

\\ ---- Test 1: Mestre_conic ground truth, I=[1,2,3,4] ----
\\ Sage docstring: Mestre_conic([1,2,3,4]) ==
\\  -2572155000*u^2 -317736000*u*v +1250755459200*v^2
\\  +2501510918400*u*w +39276887040*v*w +2736219686912*w^2
print("---- Test 1: Mestre_conic ground truth (Sage docstring, I=[1,2,3,4]) ----");
xyz1 = Mestre_xyz(1,2,3,4);
L1 = Mestre_L(xyz1[1],xyz1[2],xyz1[3]);
expected1 = [-2572155000, -317736000, 1250755459200, 2501510918400, 39276887040, 2736219686912];
got1 = [L1[1,1], 2*L1[1,2], L1[2,2], 2*L1[1,3], 2*L1[2,3], L1[3,3]];
ratios1 = vector(6, k, expected1[k]/got1[k]);
print("  x,y,z = ", xyz1);
print("  ratios expected/got (must be one constant) = ", ratios1);
print("  PASS (proportional to Sage ground truth)? ", #Set(ratios1) == 1);
print("");

\\ ---- Test 2: Mestre_conic ground truth, I=[5,6,7,8] (Sage gives x,y,z exactly) ----
print("---- Test 2: Mestre_conic ground truth (Sage docstring, I=[5,6,7,8]) ----");
xyz2 = Mestre_xyz(5,6,7,8);
print("  x,y,z = ", xyz2, "   expected = [232/1125, -1072/16875, 14695616/2109375]");
print("  PASS? ", xyz2 == [232/1125, -1072/16875, 14695616/2109375]);
L2 = Mestre_L(xyz2[1],xyz2[2],xyz2[3]);
expected2 = [-415125000, 608040000, 33065136000, 66130272000, 240829440, 10208835584];
got2 = [L2[1,1], 2*L2[1,2], L2[2,2], 2*L2[1,3], 2*L2[2,3], L2[3,3]];
ratios2 = vector(6, k, expected2[k]/got2[k]);
print("  ratios expected/got (must be one constant) = ", ratios2);
print("  PASS (proportional to Sage ground truth)? ", #Set(ratios2) == 1);
print("");

\\ ---- Test 3: full round-trip on a generic sextic over Q ----
print("---- Test 3: full round-trip, h_orig = x^6+2x^5+3x^4+5x^3+7x^2+11x+13 ----");
h_orig = x^6 + 2*x^5 + 3*x^4 + 5*x^3 + 7*x^2 + 11*x + 13;
q_orig = igusa_quadruple(h_orig);
print("  Igusa quadruple (I2,I4,I6,I10) of h_orig = ", q_orig);
f_rec = Mestre_reconstruct(q_orig[1], q_orig[2], q_orig[3], q_orig[4]);
f_rec_x = subst(f_rec, tvar, x);
q_rec = igusa_quadruple(f_rec_x);
print("  Igusa quadruple of reconstructed curve  = ", q_rec);
lam2 = q_rec[1]/q_orig[1];
ok4  = (lam2^2 == q_rec[2]/q_orig[2]);
ok6  = (lam2^3 == q_rec[3]/q_orig[3]);
ok10 = (lam2^5 == q_rec[4]/q_orig[4]);
print("  lambda^2 = ", lam2);
print("  I4 scales as lambda^4? ", ok4, "   I6 as lambda^6? ", ok6, "   I10 as lambda^10? ", ok10);
print("  ROUND-TRIP PASS (curve is Q-bar isomorphic to h_orig)? ", ok4 && ok6 && ok10);
print("");
print("  Note: reconstructed sextic has large coefficients because 'reduced'");
print("  minimisation (Sage issue #14755/#14756, still NotImplementedError");
print("  upstream) is not attempted here -- matches Sage's reduced=False path.");
print("");

\\ ============================================================
\\ Status / next step
\\ ============================================================
print("================================================================");
print("Status: Mestre Step 2 (invariants -> curve) DONE over Q.");
print("  RESEARCH_MESTRE_HOWE.md §7 gap #3 (Mestre's reconstruction,");
print("  '~500 lines; week-long effort') is closed by qfsolve+qfparam,");
print("  which PARI already provides natively for ternary forms over Q.");
print("");
print("  STILL OPEN (gap #2, the load-bearing one for secp256k1):");
print("  Igusa invariants of the ACTUAL glued surface (E x E^t)/Gamma_alpha,");
print("  as opposed to the naive product y^2=(x^3+7)(x^3+189). This script");
print("  reconstructs A CURVE from ANY given quadruple; it does not compute");
print("  the quadruple that corresponds to the Howe gluing itself.");
print("");
print("  ALSO OPEN: qfsolve needs an integer/rational matrix, so this");
print("  pipeline does not run over F_p_secp as-is. Finite-field point-");
print("  finding on the conic needs a separate diagonalize-and-take-sqrt");
print("  algorithm (isotropy of ternary forms over finite fields is");
print("  guaranteed for p odd, so a solution always exists -- just needs");
print("  an explicit constructive proof/algorithm, not Hasse-Minkowski).");
print("================================================================");

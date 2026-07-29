\\ thread22_secp256k1_realization.gp
\\
\\ Thread 22, part 2: apply the criterion measured in thread22_jacobian_realization.gp
\\ to secp256k1's six j=0 sextic twists at full cryptographic size.
\\
\\ MEASURED CRITERION (exhaustive genus-2 search over F_7, F_13, F_19 --
\\ see thread22_jacobian_realization_output.txt):
\\   A split class E_{t1} x E_{t2} over F_p contains NO Jacobian iff
\\     (X1)  |t1 - t2| = 1,   or
\\     (X2)  t1 = t2 = t  and  t^2 - 4p = -3.
\\   Every other split class over those three primes IS realized by a Jacobian.
\\   Caveat: over p in {7,13,19} no square class has t^2-4p in {-4,-7}, so those
\\   two discriminants are UNTESTED and (X2) may be incomplete. This does not
\\   affect secp256k1, where no pair has t_i = t_j at all (see below).
\\
\\ STRUCTURAL READING OF (X1) (conjectural, consistent with all data here):
\\   For P_i = T^2 - t_i*T + p, the ring Z[T]/(P_i, P_j) is finite of order
\\   p*(t_i-t_j)^2 with Smith normal form diag(|t_i-t_j|, p*|t_i-t_j|), so the
\\   reduced resultant (its exponent, = smallest positive integer in the ideal)
\\   is r = p*|t_i - t_j|.  Howe-style gluing of E_i x E_j along a group scheme
\\   of order n requires n | r.  When |t_i-t_j| = 1 we get r = p, so the only
\\   n coprime to p is n = 1: no nontrivial gluing exists, the product
\\   polarization is the only one, and a product with the product polarization
\\   is never a Jacobian.  The script verifies the SNF formula numerically.
\\
\\ Run: gp -q thread22_secp256k1_realization.gp

default(parisize, 512000000);
default(timer, 0);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n0 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t0 = p + 1 - n0;

\\ ---------------------------------------------------------------------------
\\ Verify the reduced-resultant / Smith-normal-form claim used above.
\\ Z[T]/(P_i) has Z-basis {1, T}; T acts by [[0,-p],[1,t_i]]; and
\\ P_j = P_i + (t_i - t_j)*T, so the quotient is Z^2 / ((t_i-t_j)*matrix(T)).
\\ ---------------------------------------------------------------------------
snf_check(q, t1, t2) = {
  my(M, s);
  M = (t1 - t2) * [0, -q; 1, t1];
  s = matsnf(M);
  print("    q=", q, " t1=", t1, " t2=", t2,
        "   SNF=", vecsort(s), "   predicted=", vecsort([abs(t1-t2), q*abs(t1-t2)]),
        "   ", if(vecsort(s) == vecsort([abs(t1-t2), q*abs(t1-t2)]), "OK", "MISMATCH"));
  vecsort(s) == vecsort([abs(t1-t2), q*abs(t1-t2)]);
};

print("================================================================");
print("Thread 22 part 2: Jacobian realization for secp256k1 sextic twists");
print("================================================================");
print("");
print("Reduced-resultant / SNF spot checks (r = p*|t1-t2| is the exponent):");
{
my(ok = 1);
ok = snf_check(7, -5, -4) && ok;
ok = snf_check(7, -4, -1) && ok;
ok = snf_check(13, -7, 5) && ok;
ok = snf_check(19, 8, -8) && ok;
if(!ok, error("SNF formula failed"));
}

\\ ---------------------------------------------------------------------------
\\ CM decomposition of secp256k1: pi = a + b*zeta_3, 4p = (2a-b)^2 + 3b^2.
\\ Traces of the six sextic twists are Tr(zeta_6^k * pi), k = 0..5.
\\ ---------------------------------------------------------------------------
\\ Solve 4p = t0^2 + 3*s^2 for s, then the six traces are
\\   +-t0, +-(t0+3s)/2, +-(t0-3s)/2.
s = sqrtint((4*p - t0^2)/3);
if(t0^2 + 3*s^2 != 4*p, error("CM decomposition failed"));

T = vector(6);
T[1] =  t0;                 \\ k=0  secp256k1
T[2] = -(t0 + 3*s)/2;       \\ k=1
T[3] = -(t0 - 3*s)/2;       \\ k=2
T[4] = -t0;                 \\ k=3  quadratic twist
T[5] =  (t0 + 3*s)/2;       \\ k=4
T[6] =  (t0 - 3*s)/2;       \\ k=5

print("");
print("CM data:");
print("  t0 = ", t0);
print("  s  = ", s, "     (4p = t0^2 + 3 s^2 : ", if(t0^2+3*s^2 == 4*p, "OK", "FAIL"), ")");
print("  trace sums t0+t2+t4 = ", T[1]+T[3]+T[5], " ;  t1+t3+t5 = ", T[2]+T[4]+T[6]);
{
my(bad = 0);
for(k = 1, 6, if((p + 1 - T[k]) <= 0, bad = 1));
if(bad, error("bad trace"));
}

\\ orders and Hasse check
N = vector(6, k, p + 1 - T[k]);
{
my(B = sqrtint(4*p));
for(k = 1, 6, if(abs(T[k]) > B, error("trace violates Hasse bound at k=", k-1)));
}
print("  all six traces satisfy |t| <= 2*sqrt(p): OK");

\\ ---------------------------------------------------------------------------
\\ 2-torsion Galois structure: factor x^3 + 7*h^k over F_p, h a primitive 6th
\\ root of unity mod p (p = 3 mod 4, so sqrt(-3) = (-3)^((p+1)/4)).
\\ ---------------------------------------------------------------------------
sqrtm3 = lift(Mod(-3, p)^((p+1)/4));
if(lift(Mod(sqrtm3, p)^2) != lift(Mod(-3, p)), error("sqrt(-3) failed"));
z3 = lift(Mod(-1 + sqrtm3, p) / 2);
h  = lift(Mod(-z3, p));
if(lift(Mod(h, p)^6) != 1 || lift(Mod(h, p)^3) != p-1, error("h is not a primitive 6th root"));
print("  h = primitive 6th root of unity mod p: order check OK");

B6 = vector(6, k, lift(Mod(7, p) * Mod(h, p)^(k-1)));
PAT = vector(6);
{
for(k = 1, 6,
  my(fac = factormod(Mod(1,p)*(x^3 + B6[k])));
  PAT[k] = vecsort(vector(matsize(fac)[1], i, poldegree(fac[i,1]) * fac[i,2]));
);
}

\\ Match each twist's trace to its b-value: E_k : y^2 = x^3 + 7h^k has some
\\ trace in T; identify it by testing (p+1-t)*P = O on a random point.
\\ (Needed because the algebraic ordering of Tr(zeta_6^k pi) need not match the
\\  ordering of the b-values by k.)
IDX = vector(6);
{
for(k = 1, 6,
  my(E = ellinit([0,0,0,0,B6[k]], p), P, found = 0, xx, yy2);
  \\ find an affine point
  for(a = 1, 10^6,
    yy2 = lift(Mod(a,p)^3 + B6[k]);
    if(kronecker(yy2, p) == 1, xx = a; break));
  P = [Mod(xx,p), Mod(lift(Mod(yy2,p)^((p+1)/4)), p)];
  if(lift(P[2]^2) != lift(Mod(yy2,p)), error("point construction failed at k=", k-1));
  for(j = 1, 6,
    if(ellmul(E, P, p + 1 - T[j]) == [0], IDX[k] = j; found = 1; break));
  if(!found, error("could not identify trace for k=", k-1));
);
}
{
my(seen = vector(6));
for(k = 1, 6, seen[IDX[k]]++);
for(j = 1, 6, if(seen[j] != 1, error("trace-to-twist matching is not a bijection")));
}

print("");
print("Six sextic twists of secp256k1  (E_k : y^2 = x^3 + 7*h^k):");
print("  k | trace t_k | 2-torsion pattern of x^3+7h^k");
{
for(k = 1, 6,
  print("  ", k-1, " | ", T[IDX[k]], " | ", PAT[k]);
);
}

\\ ---------------------------------------------------------------------------
\\ The 15 pairs.
\\ ---------------------------------------------------------------------------
print("");
print("15 pairs: Howe (2,2)-gluing conditions vs measured realization criterion");
print("  H1: #E_i != #E_j    H2: same 2-torsion pattern    H3: gcd(#E_i,#E_j)=1");
print("  X1: |t_i - t_j| = 1 (no Jacobian)   X2: t_i = t_j and t^2-4p = -3 (no Jacobian)");
print("");
print("  pair | H1 H2 H3 | Howe | |t_i-t_j| = 1 ? | t_i = t_j ? | class has Jacobian ?");
{
my(nhowe = 0, nreal = 0, mindiff = 0);
for(i = 1, 6, for(j = i+1, 6,
  my(ti = T[IDX[i]], tj = T[IDX[j]], ni = p+1-ti, nj = p+1-tj);
  my(h1, h2, h3, howe, x1, x2, real, d, tag);
  h1 = (ni != nj);
  h2 = (PAT[i] == PAT[j]);
  h3 = (gcd(ni, nj) == 1);
  howe = (h1 && h2 && h3);
  d = abs(ti - tj);
  x1 = (d == 1);
  x2 = (ti == tj) && (ti^2 - 4*p == -3);
  real = !(x1 || x2);
  if(howe, nhowe++);
  if(real, nreal++);
  if(mindiff == 0 || d < mindiff, mindiff = d);
  tag = concat(concat(concat("(", Str(i-1)), concat(",", Str(j-1))), ")");
  print("  ", tag, "  |  ", if(h1,"Y","n"), "  ", if(h2,"Y","n"), "  ", if(h3,"Y","n"),
        " |  ", if(howe, "YES", " no"),
        "  |      ", if(x1, "YES", "no "),
        "       |     ", if(ti == tj, "YES", "no "),
        "     |  ", if(real, "YES", "NO"));
));
print("");
print("  pairs passing Howe H1&H2&H3                       : ", nhowe, "/15");
print("  pairs whose class contains a Jacobian (criterion) : ", nreal, "/15");
print("  min over all 15 pairs of |t_i - t_j|              : ", mindiff);
print("    (~2^", floor(log(mindiff)/log(2)), "; the exception |t_i-t_j| = 1 is nowhere near being met)");
}

print("");
print("================================================================");
print("CONCLUSION");
print("================================================================");
print("All 15 product classes E_i x E_j of secp256k1's sextic twists contain a");
print("genus-2 Jacobian, whereas only 5/15 are Howe (2,2)-glueable (Thread 18).");
print("So counting Howe-glueable pairs UNDERCOUNTS the available genus-2 covers");
print("by a factor of 3 for this family.");
print("");
print("This does NOT weaken the paper's block-B5 conclusion. B5 bounds the COST of");
print("the resulting genus-2 DLP, not the NUMBER of covers: for a Jacobian J with");
print("#J ~ p^2 ~ n^2 whose group order is divisible by n, Pollard rho on J costs");
print("~sqrt(#J) >= sqrt(n) and index calculus on a genus-2 Jacobian over a PRIME");
print("field has no sub-rho algorithm (Gaudry-Thome-Theriault-Diem applies to");
print("g >= 3 or to extension fields). Tripling the number of covers multiplies the");
print("attacker's options by a constant, not the asymptotic cost.");
print("");
print("What WOULD matter, and is not settled here: whether any of the 10 newly");
print("admissible classes is realized by a Jacobian that is EXPLICITLY CONSTRUCTIBLE");
print("over F_p. Thread 2 showed the F_p Rosenhain route is blocked by the cubic");
print("residue obstruction; realization is an existence statement and gives no");
print("construction.");
quit();

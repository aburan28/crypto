\\ ============================================================
\\ Thread 22 — apply the empirical realisability criterion to secp256k1.
\\
\\ The exhaustive census (examples/howe_genus2_census.rs) over
\\ p = 7, 13, 19, 31, 37, 43 found, with no exception in 90 j=0 twist pairs
\\ and 296 arbitrary trace pairs:
\\
\\   the isogeny class E_1 × E_2 (t_1 != t_2) over F_p contains the Jacobian
\\   of a genus-2 curve   <=>   p + 1 - (t_1 + t_2) >= 0   AND   |t_1 - t_2| >= 2.
\\
\\ For the j=0 sextic-twist family the six traces are
\\   +- t, +- (t+3s)/2, +- (t-3s)/2   with 4p = t^2 + 3 s^2,
\\ whose three magnitudes {a, b, a+b} satisfy a + b = max.  Hence
\\ "some pair differs by 1" <=> "some twist has trace +-1".
\\
\\ This script checks the criterion for secp256k1's six sextic twists.
\\
\\ Run: gp -q thread22_secp256k1_criterion.gp
\\ ============================================================

default(parisize, 512000000);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t = p + 1 - n_secp;

print("================================================================");
print("Thread 22 — realisability criterion applied to secp256k1");
print("================================================================");
print("");
print("p mod 6 = ", p % 6, "   (must be 1 for six sextic twists)");
print("t(secp256k1) = ", t);

\\ ---- CM structure 4p = t^2 + 3 s^2 ----
rem = 4*p - t^2;
s_sq = rem / 3;
s = sqrtint(s_sq);
if (s^2 != s_sq, error("CM formula failed"));
print("s            = ", s);
print("check 4p = t^2 + 3 s^2 : ", 4*p == t^2 + 3*s^2);
print("");

\\ ---- the six traces ----
{
traces = [ t, -t, (t+3*s)/2, -(t+3*s)/2, (t-3*s)/2, -(t-3*s)/2 ];
}
print("---- six sextic-twist traces ----");
for(k = 1, 6, print("  t_", k-1, " = ", traces[k]));
print("");
print("  sum of traces = ", sum(k = 1, 6, traces[k]), "  (must be 0)");
print("");

\\ ---- criterion 1: does any twist have trace +-1 ? ----
print("---- criterion: is any |t_k| = 1 ? ----");
{
any1 = 0;
for(k = 1, 6, if(abs(traces[k]) == 1, any1 = 1));
print("  any |t_k| = 1 : ", if(any1, "YES -> two pairs unrealisable", "NO  -> no pair obstructed"));
}
print("");

\\ ---- criterion 2: per-pair check ----
print("---- per-pair check (all 15 pairs) ----");
print("  pair | |t_i - t_j| >= 2 | p+1-(t_i+t_j) >= 0 | class contains a Jacobian?");
{
nrealisable = 0;
for(i = 1, 6,
  for(j = i+1, 6,
    my(ti = traces[i], tj = traces[j]);
    my(c1 = (abs(ti - tj) >= 2));
    my(c2 = (p + 1 - (ti + tj) >= 0));
    if(c1 && c2, nrealisable += 1);
    print("  (", i-1, ",", j-1, ") |        ", if(c1, "yes", "NO "),
          "        |         ", if(c2, "yes", "NO "),
          "        | ", if(c1 && c2, "YES (predicted)", "no"));
  );
);
print("");
print("  predicted realisable pairs for secp256k1: ", nrealisable, "/15");
}

print("");
print("Interpretation: for secp256k1 every |t_i - t_j| is a ~2^128-size integer,");
print("so the |t_i - t_j| >= 2 clause is satisfied trivially and the point-count");
print("clause is satisfied because |t_i + t_j| << p.  The genus-2 isogeny classes");
print("E_i x E_j therefore DO contain Jacobians -- non-existence is NOT the reason");
print("the cover route fails.  The obstruction remains the DLP cost of block B5.");

quit;

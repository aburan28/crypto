\\ ============================================================
\\ Thread 19: GLV-HNP Phase 2 — actual LLL recovery on toy curve
\\ ============================================================
\\
\\ Implements the Phase 2 GLV-aware lattice attack concretely.
\\ Builds (2m+1)-dim lattice, runs qflll, checks if d is recovered.
\\
\\ Run: gp -q thread19_glv_hnp_phase2_lll.gp

default(parisize, 256000000);
default(timer, 0);

print("================================================================");
print("Thread 19: GLV-HNP Phase 2 — toy LLL recovery (n=199)");
print("================================================================");
print("");

\\ ---- Toy curve parameters ----
p_toy = 211;
b_toy = 2;
n_toy = 199;
lam = 106;    \\ GLV eigenvalue (avoid 'lambda' which is a PARI keyword)
E = ellinit([0, b_toy], p_toy);
G = ellgenerators(E)[1];
print("Toy curve: y² = x³ + ", b_toy, " / F_", p_toy);
print("n = ", n_toy, ", lam = ", lam);
print("lam² + lam + 1 mod n = ", lift(Mod(lam^2 + lam + 1, n_toy)));
print("");

\\ ---- Reproduce signatures (seed=11, matching glv_hnp_phase2_toy.gp) ----
setrand(11);
d_secret = random(n_toy - 1) + 1;
K1B = max(2, n_toy \ 100);    \\ K1_BOUND = 2 for n=199
m_sigs = 5;
print("Planted d = ", d_secret, " (UNKNOWN to attacker)");
print("K1_BOUND = ", K1B);
print("");

Avec = vector(m_sigs);
Bvec = vector(m_sigs);
k1p  = vector(m_sigs);
k2p  = vector(m_sigs);

{
cnt = 1;
while (cnt <= m_sigs,
    k1  = random(K1B - 1) + 1;
    k2  = random(n_toy - 1) + 1;
    kf  = lift(Mod(k1 + lam * k2, n_toy));
    if (kf == 0, next);
    Rpt = ellmul(E, G, kf);
    rr  = lift(Rpt[1]) % n_toy;
    if (rr == 0, next);
    hh  = random(n_toy);
    ss  = lift(Mod((hh + d_secret * rr) / kf, n_toy));
    if (ss == 0, next);
    k1p[cnt]  = k1;
    k2p[cnt]  = k2;
    Avec[cnt] = lift(Mod(hh / ss, n_toy));
    Bvec[cnt] = lift(Mod(rr / ss, n_toy));
    cnt = cnt + 1
);
}

print("Signatures (A_i, B_i) with true (k1,k2):");
{
for (ii = 1, m_sigs,
    print("  sig ", ii, ": A=", Avec[ii], " B=", Bvec[ii],
          "  [k1=", k1p[ii], " k2=", k2p[ii], "]")
);
}
print("");

\\ ---- Build (2m+1)-dim GLV lattice ----
\\
\\ Columns = lattice basis vectors (PARI convention for qflll).
\\ Lattice row index (1-based): same as column of M.
\\ We build M with BASIS VECTORS AS COLUMNS.
\\
\\ Lattice has dim = 2m+1 coordinates:
\\   coords 1..m   : k_{j,1} values (small, < K1B)
\\   coord  m+1    : d * K1B (scaled secret)
\\   coords m+2..2m+1 : k_{j,2} * K1B (scaled half-scalars)
\\
\\ Basis vectors (columns of M):
\\   Col j (j=1..m):    n * e_j   (kills j-th carry in mod-n)
\\   Col m+1 (d-col):   (B_1, ..., B_m, K1B, 0, ..., 0)
\\   Col m+1+j:         (-lam * e_j at coords 1..m; K1B at coord m+1+j)
\\
\\ Hidden vector (what we're looking for):
\\   w = (q_1, ..., q_m, d, k_{1,2}, ..., k_{m,2})
\\   where q_j = (A_j + B_j*d - lam*k_{j,2} - k_{j,1}) / n  (integer quotient)
\\
\\ M * w = (k_{1,1}, ..., k_{m,1}, d*K1B, k_{1,2}*K1B, ..., k_{m,2}*K1B)
\\ This is the "target short vector."

dim = 2 * m_sigs + 1;
M = matrix(dim, dim);

\\ Columns 1..m: n * e_j
{
for (jj = 1, m_sigs,
    M[jj, jj] = n_toy
);
}

\\ Column m+1: d-column
{
for (jj = 1, m_sigs,
    M[jj, m_sigs+1] = Bvec[jj]
);
}
M[m_sigs+1, m_sigs+1] = K1B;

\\ Columns m+2..2m+1: k_{j,2} columns
{
for (jj = 1, m_sigs,
    M[jj, m_sigs+1+jj] = -lam;
    M[m_sigs+1+jj, m_sigs+1+jj] = K1B
);
}

\\ ---- Kannan embedding: shift by -A_j in rows 1..m ----
\\
\\ We want M*w - u to be short, where u = (A_1,...,A_m, 0,...,0).
\\ Augmented matrix: embed u as extra row/column.
\\
\\ Augmented dim = dim+1:
\\   Top-left: M (basis columns)
\\   Extra column dim+1:  0 for rows 1..dim, W_embed for row dim+1
\\   Extra row dim+1: (-A_1, ..., -A_m, 0, ..., 0, W_embed)

Wemb = n_toy;    \\ weight for the embedded target
daug = dim + 1;
Maug = matrix(daug, daug);

\\ Copy M into top-left
{
for (rr = 1, dim,
    for (cc = 1, dim,
        Maug[rr, cc] = M[rr, cc]
    )
);
}

\\ Extra column (col daug): 0 for rows 1..dim
\\ (already 0 from matrix() initialisation)

\\ Extra row (row daug): -A_j for cols j=1..m_sigs; rest 0; Wemb in col daug
{
for (jj = 1, m_sigs,
    Maug[daug, jj] = -Avec[jj]
);
}
Maug[daug, daug] = Wemb;

\\ ---- Run qflll on transpose (PARI: qflll takes column vectors) ----
\\ T = qflll(Maug~) gives unimodular T with Maug~ * T LLL-reduced.
\\ Rows of reduced matrix = cols of (Maug~ * T) = rows of T~ * Maug ...
\\ Actually: reduced columns = Maug~ * T, so reduced column rr = Maug~ * T[,rr].
\\ Rows of the row-reduced matrix are rows of Maug * T~ ... let's just work
\\ with the column representation.

print("Running qflll on augmented ", daug, "x", daug, " matrix...");
T_unmod = qflll(Maug~);    \\ unimodular transform; cols of Maug~ * T are reduced
Mred_cols = Maug~ * T_unmod;    \\ each column is a reduced basis vector (coords are row indices)

print("qflll done.");
print("");

\\ Each column of Mred_cols is a lattice vector (in coords 1..daug).
\\ The last coord (coord daug) is the embedded dimension.
\\ We look for a column with abs(coord daug) = Wemb.

print("Scanning reduced basis columns for embedded target (coord_", daug, " = ±", Wemb, ")...");
print("");
found_d = -1;
{
for (cc = 1, daug,
    ecol = Mred_cols[, cc];    \\ column cc as a vector
    ec_last = ecol[daug];
    if (abs(ec_last) != Wemb, next);
    \\ This column has the embedding marker
    sign_v = if(ec_last > 0, 1, -1);
    k1_vals = vector(m_sigs, jj, sign_v * ecol[jj]);
    d_raw   = sign_v * ecol[m_sigs+1];
    k2_vals = vector(m_sigs, jj, sign_v * ecol[m_sigs+1+jj]);
    print("Found candidate column ", cc, ":");
    print("  k1 values:  ", k1_vals);
    print("  d_raw (= d*K1B): ", d_raw, "  K1B = ", K1B);
    print("  k2 values:  ", k2_vals);
    \\ Check k1 values in [0, K1B)
    k1_ok = 1;
    for (jj = 1, m_sigs,
        if (k1_vals[jj] < 0 || k1_vals[jj] >= K1B, k1_ok = 0)
    );
    if (!k1_ok,
        print("  k1 values out of [0, K1B) — not the target.");
        next
    );
    if (d_raw % K1B != 0,
        print("  d_raw not divisible by K1B — not the target.");
        next
    );
    d_cand = d_raw / K1B;
    d_cand = lift(Mod(d_cand, n_toy));
    print("  d_candidate = ", d_cand, "  planted d = ", d_secret);
    if (d_cand == d_secret,
        print("  *** d RECOVERED! *** d = ", d_cand, " ✓");
        found_d = d_cand
    ,
        print("  Candidate d ≠ planted d (checking all cols...)")
    )
);
}

if (found_d == -1,
    print("Direct column scan: d not found.  Printing first m_sigs+3 coords of all cols:");
    for (cc = 1, daug,
        ecol = Mred_cols[, cc];
        print("  col ", cc, ": last=", ecol[daug], " d_pos=", ecol[m_sigs+1],
              " k1s=[", ecol[1], ",", ecol[2], ",", ecol[3], ",", ecol[4], ",", ecol[5], "]")
    )
);

print("");
print("================================================================");
print("Thread 19 LLL attempt complete.");
print("================================================================");

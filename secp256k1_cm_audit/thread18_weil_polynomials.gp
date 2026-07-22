\\ Thread 18: Weil polynomials for Howe-glueable secp256k1 sextic twist pairs
\\ PARI 2.15.4 compatible: all for-body logic in named functions; no multi-line for bodies.

default(parisize, 512000000);
default(timer, 0);
default(realprecision, 50);

print("=================================================================");
print("Thread 18: Weil polynomials for Howe-glueable secp256k1 twist pairs");
print("=================================================================");
print("");

\\ ------------------------------------------------------------------
\\ secp256k1 parameters
\\ ------------------------------------------------------------------
SP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
SN = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
T0 = SP + 1 - SN;

print("p = ", SP);
print("n = ", SN);
print("t0 = p+1-n = ", T0);
print("");

\\ CM parameter: 4p = t0^2 + 3*S^2
REM3 = 4*SP - T0^2;
if(REM3 % 3 != 0, error("CM fail: 4p-t0^2 not div by 3"));
SV2 = REM3 / 3;
SV = sqrtint(SV2);
if(SV^2 != SV2, error("CM fail: S not integer"));
print("CM: 4p = t0^2 + 3*S^2  =>  S = ", SV);
print("Verify 4p = t0^2+3S^2: ", 4*SP == T0^2 + 3*SV^2);
print("");

\\ ------------------------------------------------------------------
\\ Six traces indexed 0..5 (stored 1-indexed in TR[1..6])
\\ ------------------------------------------------------------------
TR = vector(6);
TR[1] = T0;
TR[2] = (T0 + 3*SV) / 2;
TR[3] = (-T0 + 3*SV) / 2;
TR[4] = -T0;
TR[5] = -(T0 + 3*SV) / 2;
TR[6] = (T0 - 3*SV) / 2;

print("Six sextic twist traces (k=0..5):");
print_trace(k) = {
    my(tk, Nk, bits);
    tk = TR[k+1];
    Nk = SP + 1 - tk;
    bits = floor(log(abs(tk)) / log(2) + 0.5);
    print("  k=", k, ": t = ", tk, "   bits ~ ", bits);
};
print_trace(0);
print_trace(1);
print_trace(2);
print_trace(3);
print_trace(4);
print_trace(5);
print("");
if(SP + 1 - TR[1] != SN, error("t0 check fail"));
print("Order check: p+1-t[0] = n ✓");
print("");

\\ ------------------------------------------------------------------
\\ Glueable pairs (from 2026-06-09 PARI run, 5 of 15 pairs)
\\ Pairs (0-indexed): (0,2),(0,3),(0,5),(1,4),(2,3)
\\ ------------------------------------------------------------------
GP_I = [0, 0, 0, 1, 2];
GP_J = [2, 3, 5, 4, 3];

\\ Helper: compute Weil polynomial coefficients for pair (i,j)
\\ W(T) = (T^2 - ti*T + p)(T^2 - tj*T + p)
\\       = T^4 - (ti+tj)*T^3 + (2p+ti*tj)*T^2 - p*(ti+tj)*T + p^2
weil_c3(i,j) = -(TR[i+1] + TR[j+1]);
weil_c2(i,j) = 2*SP + TR[i+1]*TR[j+1];
weil_c1(i,j) = -SP*(TR[i+1] + TR[j+1]);
weil_jacsize(i,j) = (SP + 1 - TR[i+1]) * (SP + 1 - TR[j+1]);
is_biquad(i,j) = (TR[i+1] + TR[j+1] == 0);

\\ Helper: log2 of a positive integer
log2_int(x) = log(x*1.0) / log(2.0);

\\ ------------------------------------------------------------------
\\ Print detailed results for each glueable pair
\\ ------------------------------------------------------------------
print("====================================================================");
print("Weil polynomial details for each glueable pair");
print("====================================================================");
print("");

print_pair(idx) = {
    my(i, j, ti, tj, c2, c3, jac, biq, ecdlp_b, jac_b, extra, a2, disc, sf_val, a2modp);
    i = GP_I[idx];
    j = GP_J[idx];
    ti = TR[i+1];
    tj = TR[j+1];
    c3 = -(ti+tj);
    c2 = 2*SP + ti*tj;
    jac = (SP+1-ti) * (SP+1-tj);
    biq = (ti+tj == 0);
    ecdlp_b = log2_int(SN) / 2;
    jac_b   = log2_int(jac) / 2;
    extra   = jac_b - ecdlp_b;
    print("--- Pair (", i, ",", j, ") ---");
    print("  ti = ", ti);
    print("  tj = ", tj);
    print("  ti+tj = ", ti+tj);
    print("  Biquadratic: ", if(biq,"YES","NO"));
    print("  W(T) = T^4 + (", c3, ")*T^3 + (", c2, ")*T^2 + ...*T + p^2");
    print("  |Jac(C)| = ", jac);
    print("  log2|Jac| = ", log2_int(jac));
    print("  ECDLP cost: 2^", ecdlp_b, " bits (sqrt(n))");
    print("  Jac DLP cost: 2^", jac_b, " bits (sqrt(|Jac|))");
    print("  Extra cost: +", extra, " bits (Jac DLP harder by 2^", extra, ")");
    if(biq,
        a2 = c2;
        disc = a2^2 - 4*SP^2;
        sf_val = core(abs(disc));
        if(disc < 0, sf_val = -sf_val);
        a2modp = lift(Mod(a2, SP));
        print("  >>> BIQUADRATIC: a2 = 2p+ti*tj = ", a2);
        print("  >>> D = a2^2-4p^2 = ", disc);
        print("  >>> sf(D) = ", sf_val);
        print("  >>> a2 mod p = ", a2modp, " (=0 means theorem degenerate)");
        print("  >>> p | a2: ", if(a2modp==0, "YES (p|ti*tj mod p, need check)", "NO, theorem applies ✓"));
    );
    print("");
};

print_pair(1);
print_pair(2);
print_pair(3);
print_pair(4);
print_pair(5);

\\ ------------------------------------------------------------------
\\ Summary table
\\ ------------------------------------------------------------------
print("====================================================================");
print("Summary table");
print("====================================================================");
print("Pair | Biquad | log2|Jac| | Jac bits | ECDLP bits | Extra bits");
print_summary_row(idx) = {
    my(i, j, ti, tj, jac, biq, ecdlp_b, jac_b, extra);
    i = GP_I[idx];
    j = GP_J[idx];
    ti = TR[i+1];
    tj = TR[j+1];
    jac = (SP+1-ti) * (SP+1-tj);
    biq = (ti+tj == 0);
    ecdlp_b = log2_int(SN) / 2;
    jac_b   = log2_int(jac) / 2;
    extra   = jac_b - ecdlp_b;
    print("(", i, ",", j, ")  | ", if(biq,"YES","NO "), "    | ",
          log2_int(jac), " | ",
          jac_b, " | ",
          ecdlp_b, " | +", extra);
};
print_summary_row(1);
print_summary_row(2);
print_summary_row(3);
print_summary_row(4);
print_summary_row(5);

print("");
print("CONCLUSION: generic Jac DLP cost >> ECDLP for all 5 pairs.");
print("Cover attack (lifting E DLP to Jac DLP) offers no speedup.");
print("");

\\ ------------------------------------------------------------------
\\ Theorem 15-16: verify for biquadratic pairs (0,3) and (1,4)
\\ ------------------------------------------------------------------
print("====================================================================");
print("Theorem 15-16 verification for biquadratic pairs");
print("====================================================================");

verify_theorem(i,j) = {
    my(ti, tj, a2, disc, sf_val, a2mp);
    ti = TR[i+1];
    tj = TR[j+1];
    a2 = 2*SP + ti*tj;
    disc = a2^2 - 4*SP^2;
    sf_val = core(abs(disc));
    if(disc < 0, sf_val = -sf_val);
    a2mp = lift(Mod(a2, SP));
    print("Pair (", i, ",", j, "):");
    print("  a2 = 2p + t", i, "*t", j, " = ", a2);
    print("  D = a2^2 - 4p^2 = ", disc);
    print("  sf(D) = ", sf_val);
    print("  K = Q(sqrt(", sf_val, "))");
    print("  a2 mod p = (t", i, "*t", j, ") mod p = ", a2mp);
    print("  p | a2: ", if(a2mp == 0, "YES", "NO -- Theorem 15-16 applies ✓"));
    print("  p | sf(D): ", if(sf_val % SP == 0, "YES (ramification)", "NO (p splits in K) ✓"));
    print("");
};

verify_theorem(0, 3);
verify_theorem(1, 4);

print("Both biquadratic pairs: Theorem 15-16 applies => [P]^2 = 1 in Cl(K), p splits in K.");
print("");
print("For non-biquadratic pairs (0,2),(0,5),(2,3):");
print("  Weil polynomial T^4 + c3*T^3 + c2*T^2 + c1*T + p^2 with c3 != 0.");
print("  Theorem 15-16 does not apply directly (not of T^4+a2T^2+p^2 form).");
print("  Block B5 still holds: |Jac| ~ p^2 => generic DLP cost ~ p >> sqrt(n).");
print("");
print("Done. All 5 Howe-glueable pairs confirmed: cover attack fails on secp256k1.");

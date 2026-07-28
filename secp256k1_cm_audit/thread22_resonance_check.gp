\\ thread22_resonance_check.gp
\\ Validate resonance obstruction over multiple proxy primes
\\ Hypothesis: pair (i,j) is NEGATIVE in Z/3Z search iff t_i+t_j equals some trace t_k
default(parisize, 256000000); default(timer, 0);

check_prime(p_in) = {
    local(t0, s0, T_v, N_v, u_val, b_v, dp, tr_set);
    local(qi, qj, ti, tj, sm, nm, a2, a1, tcp_ok);

    \\ Find j=0 curve with prime order at this prime if possible, else just use y^2=x^3+b_try
    \\ For simplicity: find b s.t. y^2=x^3+b has prime-ish order and CM
    \\ Use b=7 always (the secp256k1 b-value), trust CM formula works

    \\ Detect the trace of E: y^2=x^3+7 over F_p
    local(E0, n00, t00, sq_test);
    E0 = ellinit([0, 7], p_in);
    n00 = ellcard(E0);
    t00 = p_in + 1 - n00;
    \\ CM: 4p = t^2 + 3s^2
    local(rem, sq, s_try);
    rem = 4*p_in - t00^2;
    if (rem % 3 != 0 || rem < 0, print("  p=", p_in, ": 4p-t^2 not divisible by 3, skip"); return(0));
    sq = rem / 3;
    s_try = sqrtint(sq);
    if (s_try^2 != sq, print("  p=", p_in, ": s not integer, skip"); return(0));

    t0 = t00; s0 = s_try;
    T_v = [t0, (t0-3*s0)/2, -(t0+3*s0)/2, -t0, (3*s0-t0)/2, (t0+3*s0)/2];
    \\ Check all are integers
    if (denominator(T_v[2]) != 1, print("  p=", p_in, ": odd trace halving, skip"); return(0));

    N_v = vector(6, k, p_in + 1 - T_v[k]);

    \\ 6th root of unity
    u_val = lift(znprimroot(p_in) ^ ((p_in - 1) / 6));
    b_v = vector(6, k, lift(Mod(7, p_in) * Mod(u_val, p_in)^(k-1)));

    \\ 2-torsion patterns
    dp = vector(6);
    for (kk = 1, 6,
        local(f1, fac1, nr1, pat1);
        f1 = Mod(1, p_in) * (x^3 + b_v[kk]);
        fac1 = factor(f1);
        nr1 = matsize(fac1)[1];
        pat1 = [];
        for (i = 1, nr1, pat1 = concat(pat1, [poldegree(fac1[i,1])]));
        dp[kk] = pat1
    );

    \\ Trace set
    tr_set = Set(T_v);

    print("p=", p_in, "  traces=", T_v, "  s=", s0);

    \\ Find qualifying pairs and check each
    for (ii = 1, 6,
        for (jj = ii+1, 6,
            local(h1, h2, h3);
            h1 = (N_v[ii] != N_v[jj]);
            h2 = (dp[ii] == dp[jj]);
            h3 = (gcd(N_v[ii], N_v[jj]) == 1);
            if (!(h1 && h2 && h3), next);
            \\ Qualifying pair
            ti = T_v[ii]; tj = T_v[jj];
            \\ Resonance check: is ti+tj in the trace set?
            local(resonance_flag);
            resonance_flag = setsearch(tr_set, ti+tj) > 0;
            \\ Predicted outcome
            local(predicted);
            predicted = if(resonance_flag, "NEGATIVE", "POSITIVE");
            \\ Actual Z/3Z search
            a2 = ti*tj + 2*p_in;
            a1 = -p_in*(ti+tj);
            nm = 0; sm = 0;
            for (a = 0, p_in-1,
                for (b = 0, p_in-1,
                    local(f2, dm, cp2, c3, c2_v, c1, c0);
                    f2 = x^6 + a*x^3 + b;
                    dm = lift(poldisc(f2) * Mod(1, p_in));
                    if (dm == 0, next);
                    sm++;
                    cp2 = hyperellcharpoly(f2 * Mod(1, p_in));
                    c3 = polcoef(cp2, 3); c2_v = polcoef(cp2, 2);
                    c1 = polcoef(cp2, 1); c0 = polcoef(cp2, 0);
                    if (c3 == -(ti+tj) && c2_v == a2 && c1 == a1 && c0 == p_in^2,
                        nm++
                    )
                )
            );
            local(actual);
            actual = if(nm > 0, "POSITIVE", "NEGATIVE");
            local(match_flag);
            match_flag = if(actual == predicted, "OK", "MISMATCH!");
            print("  pair(", ii-1, ",", jj-1, ") t=(", ti, ",", tj, ") sum=", ti+tj,
                  " resonance=", resonance_flag, " predicted=", predicted,
                  " actual=", actual, "(", nm, "/", sm, ") ", match_flag)
        )
    );
    print("")
};

print("=================================================================");
print("Resonance hypothesis validation across proxy primes");
print("H: Z/3Z search NEGATIVE iff t_i+t_j is a trace of another twist");
print("=================================================================");
print("");
check_prime(43);
check_prime(67);
check_prime(97);
print("=================================================================");
print("Done.");

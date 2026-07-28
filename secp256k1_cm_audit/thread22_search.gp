\\ thread22_search.gp -- direct Z/3Z search for qualifying pairs over p=43
default(parisize, 256000000); default(timer, 0);
p = 43;

\\ 5 qualifying pairs and their target char polys (from thread22_z3z analysis):
\\ Traces over p=43: T=[13,5,-8,-13,-5,8], same-dp pairs:
\\ (0,2) t=(13,-8)  (0,3) t=(13,-13)  (0,5) t=(13,8)  (1,4) t=(5,-5)  (2,3) t=(-8,-13)

search_pair(ti, tj, label) = {
    local(a2, a1, a0, n_match, n_sm);
    \\ target: (T^2-ti*T+p)(T^2-tj*T+p)
    a2 = ti*tj + 2*p;
    a1 = -p*(ti+tj);
    a0 = p^2;
    print("--- Pair ", label, ": t=(", ti, ",", tj, ") ---");
    print("  target charpoly = x^4 + ", -(ti+tj), "*x^3 + ", a2, "*x^2 + ", a1, "*x + ", a0);
    n_match = 0; n_sm = 0;
    for (a = 0, p-1,
        for (b = 0, p-1,
            local(f, dm, cp, c3, c2, c1, c0, c4);
            f = x^6 + a*x^3 + b;
            dm = lift(poldisc(f) * Mod(1, p));
            if (dm == 0, next);
            n_sm++;
            cp = hyperellcharpoly(f * Mod(1, p));
            \\ cp coefficients (degree 4): cp = x^4 + c3*x^3 + c2*x^2 + c1*x + c0
            c3 = polcoef(cp, 3); c2 = polcoef(cp, 2); c1 = polcoef(cp, 1); c0 = polcoef(cp, 0);
            if (c3 == -(ti+tj) && c2 == a2 && c1 == a1 && c0 == a0,
                n_match++;
                if (n_match <= 3, print("  match a=", a, " b=", b, " cp=", cp))
            )
        )
    );
    print("  smooth=", n_sm, "  matches=", n_match);
    if (n_match > 0, print("  POSITIVE"), print("  NEGATIVE"));
    print("")
};

print("Thread 22 -- Z/3Z search for qualifying pairs, p=43");
print("================================================================");
print("");

search_pair(13, -13, "(0,3) quad-twist");
search_pair(13, -8,  "(0,2) non-quad");
search_pair(13, 8,   "(0,5) non-quad");
search_pair(5, -5,   "(1,4) quad-twist TypeII");
search_pair(-8, -13, "(2,3) non-quad");

print("================================================================");
print("Done.");

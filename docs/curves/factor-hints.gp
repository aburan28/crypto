\\ Regenerate docs/curves/factor-hints.txt (see docs/curves/TRAITS.md, "Factor hints").
\\
\\   ./target/release/curve_traits unfactored > /tmp/unfactored.txt
\\   UNFACTORED=/tmp/unfactored.txt HINT_SECONDS=60 \
\\     gp -q -s 1000000000 docs/curves/factor-hints.gp < /dev/null > /tmp/hints.txt
\\   mv /tmp/hints.txt docs/curves/factor-hints.txt
\\
\\ Each composite gets HINT_SECONDS of PARI's factor(); one that does not finish
\\ contributes nothing.  The output keeps every existing hint and adds each new
\\ prime factor above 2^16, ascending, one per line.  curve_traits re-checks
\\ every hint (Miller-Rabin) and only ever uses one as a trial divisor, so a
\\ hint can speed a build up but cannot change what is proved.
U = readvec(getenv("UNFACTORED"));
H = iferr(readvec("docs/curves/factor-hints.txt"), e, []);
T = eval(getenv("HINT_SECONDS"));
P = Set(H);
{
  for (i = 1, #U,
    \\ alarm() returns an error object, not a matrix, on timeout.
    f = iferr(alarm(T, factor(U[i])), e, 0);
    if (type(f) == "t_MAT",
      P = setunion(P, Set(select(p -> p > 2^16 && ispseudoprime(p), f[, 1]~)))));
}
for (i = 1, #P, print(P[i]));
quit;

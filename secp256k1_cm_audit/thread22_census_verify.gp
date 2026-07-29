\\ ============================================================
\\ Thread 22 — independent PARI cross-check of the Rust genus-2 census
\\   examples/howe_genus2_census.rs
\\
\\ Two independent checks:
\\  (A) For each curve the Rust census reported as realising P_i·P_j,
\\      recompute charpoly(Frob | Jac C) with PARI's hyperellcharpoly
\\      and compare against (T²-t_i T+p)(T²-t_j T+p).
\\  (B) For p = 7 re-derive the FULL set of realisable genus-2 Frobenius
\\      char polys by brute force over an independent model family
\\      (all squarefree sextics a6 x^6 + ... + a0, a6 != 0, plus all
\\      squarefree quintics), and confirm the two pairs the Rust census
\\      reports as unrealisable are genuinely absent.
\\
\\ Run: gp -q thread22_census_verify.gp
\\ ============================================================

default(parisize, 512000000);

\\ ---- sextic twists of y^2 = x^3 + b over F_p, p = 1 mod 6 ------------
twist_traces(p) =
{
  my(g, tw, b, E);
  g = znprimroot(p);
  tw = vector(6);
  for(k = 0, 5,
    b = lift(g^k);
    E = ellinit([0, b], p);
    tw[k+1] = ellap(E);
  );
  tw;
}

print("================================================================");
print("Thread 22 verification — PARI cross-check");
print("================================================================");

\\ ---- (A) verify the census hits -------------------------------------
\\ Format: [p, i, j, f(x)]  (i, j are 0-based twist indices)
\\ NOTE: PARI 2.15 parses a script line-by-line, so a multi-line `[` ... `]`
\\ literal must be wrapped in braces (same gotcha as howe_22_kernel.gp).
{
hits = [
  [7,  1, 2, x^6 + 3],
  [7,  1, 4, x^6 + 2*x^3 + 6],
  [7,  1, 5, 3*x^6 + 6],
  [7,  2, 5, x^6 + 2*x^3 + x + 6],
  [7,  4, 5, 3*x^6 + 2],
  [7,  0, 3, x^6 + 6],
  [7,  0, 2, x^6 + 2*x^3 + x^2 + 4],
  [13, 1, 2, x^6 + 2],
  [13, 1, 4, x^6 + x^2 + 5],
  [13, 1, 5, 2*x^6 + 6],
  [13, 4, 5, 2*x^6 + 4],
  [13, 2, 5, x^6 + 4*x^3 + 5],
  [13, 0, 5, x^6 + x^3 + 11]
];
}

print("");
print("---- (A) recompute charpoly(Jac C) for reported census hits ----");
print("");
print("   p | (i,j) | f(x)                       | PARI charpoly(Jac)          | match");
nA = 0; okA = 0;
{
for(h = 1, #hits,
  my(p = hits[h][1], i = hits[h][2], j = hits[h][3], f = hits[h][4]);
  my(tw = twist_traces(p));
  my(ti = tw[i+1], tj = tw[j+1]);
  my(target = (x^2 - ti*x + p) * (x^2 - tj*x + p));
  my(got = hyperellcharpoly(Mod(1,p)*f));
  my(ok = (got == target));
  nA += 1; if(ok, okA += 1);
  print("  ", p, " | (", i, ",", j, ") | ", f, "  |  ", got, "  |  ", if(ok, "OK", "*** MISMATCH ***"));
);
}
print("");
print("  (A) ", okA, "/", nA, " hits confirmed by PARI hyperellcharpoly");

\\ ---- (B) independent exhaustive p = 7 census ------------------------
print("");
print("---- (B) independent exhaustive census over F_7 ----");
print("");
{
  my(p = 7, seen = List(), tw, s1, s2, target, f, cp, realised);
  tw = twist_traces(p);
  print("  twist traces over F_7: ", tw);

  \\ all sextics a6*x^6 + ... + a0 with a6 != 0, squarefree
  for(a6 = 1, p-1,
   for(a5 = 0, p-1,
    for(a4 = 0, p-1,
     for(a3 = 0, p-1,
      for(a2 = 0, p-1,
       for(a1 = 0, p-1,
        for(a0 = 0, p-1,
          f = Mod(1,p)*(a6*x^6 + a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0);
          if(poldisc(f) != Mod(0,p),
            cp = hyperellcharpoly(f);
            listput(seen, cp);
          );
        )))))));
  \\ all quintics a5*x^5 + ... + a0 with a5 != 0, squarefree
  for(a5 = 1, p-1,
   for(a4 = 0, p-1,
    for(a3 = 0, p-1,
     for(a2 = 0, p-1,
      for(a1 = 0, p-1,
       for(a0 = 0, p-1,
         f = Mod(1,p)*(a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0);
         if(poldisc(f) != Mod(0,p),
           cp = hyperellcharpoly(f);
           listput(seen, cp);
         );
       ))))));

  seen = Set(Vec(seen));
  print("  distinct genus-2 Frobenius char polys over F_7: ", #seen);
  print("");
  print("  pair | s1 | s2  | realised (PARI exhaustive)");
  for(i = 0, 5,
    for(j = i+1, 5,
      my(ti = tw[i+1], tj = tw[j+1]);
      target = (x^2 - ti*x + p) * (x^2 - tj*x + p);
      s1 = ti + tj;
      s2 = ti*tj + 2*p;
      realised = setsearch(seen, target);
      print("  (", i, ",", j, ") | ", s1, " | ", s2, "  | ", if(realised, "YES", "NO"));
    );
  );
}

print("");
print("Done.");
quit;

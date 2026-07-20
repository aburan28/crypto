\\ Thread 17: Verify p always splits (not inert, not ramified) in K=Q(sqrt(sf))
\\ when p is odd prime and p does not divide a2.
\\
\\ Algebraic argument:
\\   D = a2^2 - 4*p^2,  sf = squarefree(D).
\\   Inert impossible: Thread 16 showed the only norm-p^2 ideal when p inert is (p),
\\     but (beta)=(p) forces p|a2, excluded.
\\   Ramified impossible: p ramifies in Q(sqrt(sf)) iff p|disc(K).
\\     For odd p: disc(K) = sf or 4*sf, so p|disc(K) iff p|sf.
\\     sf divides D = a2^2 - 4p^2 ≡ a2^2 (mod p).
\\     If p|sf then p|D hence p|a2 -- excluded.
\\   Therefore p SPLITS.

check1(p, a2) =
{
  local(D, sf, discK, nf, fact, bnf, idealP, cls, cyc, expv, ord2ok, i);

  D = a2^2 - 4*p^2;
  sf = core(D);

  \\ discriminant of Q(sqrt(sf))
  if (sf % 4 == 1, discK = sf, discK = 4*sf);

  \\ Ramification check: p should NOT divide discK
  if (discK % p == 0,
    print("  RAM FAIL at p=",p," a2=",a2);
    return(0)
  );

  \\ Build NF and factor p
  nf = nfinit(x^2 - sf);
  fact = idealprimedec(nf, p);
  if (#fact != 2,
    print("  SPLIT FAIL at p=",p," a2=",a2," #primes=",#fact);
    return(0)
  );

  \\ Check [P]^2 = 1 in Cl(K)
  bnf = bnfinit(x^2 - sf, 1);
  idealP = bnfisprincipal(bnf, fact[1], 1);
  cls = bnf.clgp;
  cyc = cls.cyc;
  expv = idealP[1];
  ord2ok = 1;
  for (i = 1, #cyc,
    if ((2 * expv[i]) % cyc[i] != 0, ord2ok = 0)
  );
  if (!ord2ok,
    print("  ORDER2 FAIL at p=",p," a2=",a2);
    return(0)
  );

  print("  p=",p,"  a2=",a2,"  sf=",sf,"  disc=",discK,"  h=",cls.no,"  PASS");
  return(1);
}

print("=== Thread 17: p always splits, [P]^2=1 (10 cases) ===");
print("");

total = 0; passed = 0;

\\ Case 1: p=7, a2=3
total++; passed += check1(7, 3);
\\ Case 2: p=13, a2=5
total++; passed += check1(13, 5);
\\ Case 3: p=19, a2=7
total++; passed += check1(19, 7);
\\ Case 4: p=31, a2=11
total++; passed += check1(31, 11);
\\ Case 5: p=41, a2=13
total++; passed += check1(41, 13);
\\ Case 6: p=47, a2=22  (Thread 16 case, sf=-58, m=12)
total++; passed += check1(47, 22);
\\ Case 7: p=53, a2=10  (Thread 16 case, h=12)
total++; passed += check1(53, 10);
\\ Case 8: p=61, a2=17
total++; passed += check1(61, 17);
\\ Case 9: p=71, a2=8
total++; passed += check1(71, 8);
\\ Case 10: p=89, a2=21
total++; passed += check1(89, 21);

print("");
print(passed,"/",total," cases passed.");
if (passed == total,
  print("ALL PASS: p splits (not inert, not ramified); [P]^2=1 -- algebraic proof confirmed");
  ,
  print("FAILURE: see cases above")
);

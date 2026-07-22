\\ thread18_splitting_split_confirm.gp
\\ Verify that p ALWAYS SPLITS (never ramifies, never is inert) in Q(sqrt(sf))
\\ when p is odd prime and p∤a₂, where sf = core(a₂²-4p²).
\\
\\ Algebraic proof (from Thread 17, Proposition prop:biquadratic-order2):
\\   x²+a₂x+p² ≡ x·(x+a₂) mod p  (since p²≡0)
\\   roots {0, -a₂} are distinct iff a₂≢0 (mod p) ↔ p∤a₂  ✓
\\   Two distinct roots ↔ polynomial splits mod p ↔ p splits in Q[x]/(x²-sf)
\\
\\ This script verifies numerically over 20 (p,a₂) pairs (mixed D<0 and D>0).

check_case(p, a2) =
{
  local(D, sf, m, kron, root0, rootA, ok_noram, ok_split, ok_roots);

  \\ Precondition
  if(a2 % p == 0, error("p|a2: excluded by hypothesis"));

  D  = a2^2 - 4*p^2;
  sf = core(D);
  m  = sqrtint(D / sf);   \\ D = m^2 * sf

  \\ (1) No ramification: p must NOT divide sf.
  \\ Since v_p(D)=0 (because D≡a₂²≢0 mod p) and sf|D with sf squarefree,
  \\ any prime dividing sf also divides D, but p∤D so p∤sf.
  ok_noram = (sf % p != 0);

  \\ (2) Splitting: Kronecker symbol (sf/p) must equal 1.
  kron = kronecker(sf, p);
  ok_split = (kron == 1);

  \\ (3) Direct root check: x²+a₂x+p² ≡ 0 mod p at x=0 and x=-a₂.
  root0 = Mod(0,p)^2 + Mod(a2,p)*Mod(0,p) + Mod(p^2,p);
  rootA = Mod(-a2,p)^2 + Mod(a2,p)*Mod(-a2,p) + Mod(p^2,p);
  ok_roots = (root0 == 0) && (rootA == 0);

  return([ok_noram, ok_split, ok_roots, sf, kron]);
}

{
  \\ 20 test cases: 10 with D<0 (imaginary quadratic K), 10 with D>0 (real K)
  \\ Format: [p, a2] where p prime, p∤a2, and D = a2^2-4p^2

  \\ D < 0 pairs (a2 < 2p)
  pairs_neg = [[7,3],[11,4],[13,5],[17,6],[19,8],[23,9],[29,12],[31,7],[37,14],[41,16]];

  \\ D > 0 pairs (a2 > 2p)
  pairs_pos = [[5,11],[7,15],[11,23],[13,27],[17,35],[19,39],[23,47],[29,59],[31,63],[37,75]];

  all_pairs = concat(pairs_neg, pairs_pos);

  pass = 0; fail = 0;

  print("=== Thread 18: Splitting verification (p∤a₂ → p splits in Q(√sf)) ===");
  print("");
  printf("%-4s  %-5s  %-8s  %-6s  %-7s  %-7s  %-7s\n",
         "p", "a2", "sf", "kron", "no-ram", "split", "roots");
  print("---------------------------------------------------------------");

  for(i = 1, #all_pairs,
    p  = all_pairs[i][1];
    a2 = all_pairs[i][2];

    res = check_case(p, a2);
    ok_noram = res[1];
    ok_split = res[2];
    ok_roots = res[3];
    sf       = res[4];
    kron     = res[5];

    ok = ok_noram && ok_split && ok_roots;
    if(ok, pass++, fail++);

    printf("%-4d  %-5d  %-8d  %-6d  %-7s  %-7s  %-7s  %s\n",
           p, a2, sf, kron,
           if(ok_noram,"OK","FAIL"),
           if(ok_split,"OK","FAIL"),
           if(ok_roots,"OK","FAIL"),
           if(ok,"✓","✗"));
  );

  print("---------------------------------------------------------------");
  printf("Result: %d/%d passed\n", pass, pass+fail);
  if(fail == 0,
    print("PASS: p ALWAYS SPLITS in Q(sqrt(sf)) when p∤a₂  [Thread 17 Prop confirmed]"),
    print("FAIL: ",fail," cases did not pass — check above")
  );

  \\ -------------------------------------------------------------------------
  \\ Extra: show that a ramified case WOULD require p | a2
  \\ (i.e., if a2 = p*k for k∈Z, then sf can be divisible by p)
  print("");
  print("=== Ramification sanity check: p|a₂ case ===");
  p_ex = 7; a2_ex = 7*2;  \\ p|a2: a2=14, p=7
  D_ex = a2_ex^2 - 4*p_ex^2;
  sf_ex = core(D_ex);
  printf("p=%d a2=%d (p|a2): D=%d sf=%d  p|sf? %s\n",
         p_ex, a2_ex, D_ex, sf_ex, if(sf_ex % p_ex == 0, "YES (ramified!)", "no"));

  p_ex2 = 11; a2_ex2 = 11*3;
  D_ex2 = a2_ex2^2 - 4*p_ex2^2;
  sf_ex2 = core(D_ex2);
  printf("p=%d a2=%d (p|a2): D=%d sf=%d  p|sf? %s\n",
         p_ex2, a2_ex2, D_ex2, sf_ex2, if(sf_ex2 % p_ex2 == 0, "YES (ramified!)", "no"));
}

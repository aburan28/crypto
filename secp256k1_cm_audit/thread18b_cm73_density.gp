\\ thread18b_cm73_density.gp
\\ NEW ATTACK VARIANT: CM-73 density search
\\
\\ Prior work (Threads 11-14) showed that sf(a₂²-4p²) = -219 ONLY for primes
\\ p in "norm form" p = k²+219 with a₂ = -2k (within k≤199). Those 4 primes
\\ {19,37,79,109} were the only norm-form examples found.
\\
\\ Research question: for ARBITRARY primes p, how often does there exist
\\ an a₂ with core(a₂²-4p²) = -219? If such pairs are widespread, it suggests
\\ the CM-73 cover-attack window is wider than previously known.
\\
\\ Method: a₂²+219m² = 4p². Iterate m from 1 to floor(2p/√219), check
\\ if 4p²-219m² is a perfect square. Then verify core.
\\
\\ For p in [100, 5000]: full search. For p in [10000, 20000]: sampled run.

print("=== Thread 18b: CM-73 density search (sf = -219) ===");
print("");

{
  target_sf = -219;
  sqr219 = sqrt(219);
  count_small = 0;
  norm_form_hits = List();
  non_norm_hits  = List();

  \\ Full search over p in [5, 5000]
  forprime(p = 5, 5000,
    mmax = floor(2*p / sqr219);
    for(m = 1, mmax,
      rhs = 4*p^2 - 219*m^2;
      if(rhs <= 0, break);
      a2sq = rhs;
      a2 = sqrtint(a2sq);
      if(a2^2 != a2sq, next);
      \\ We have a2^2 + 219*m^2 = 4p^2; verify sf
      D = a2^2 - 4*p^2;
      sf = core(D);
      if(sf == target_sf,
        count_small++;
        \\ Check norm-form: p = k^2 + 219 for some k
        k = sqrtint(p - 219);
        is_nf = (k^2 + 219 == p);
        if(is_nf,
          listput(norm_form_hits, [p, a2, m, k]),
          listput(non_norm_hits,  [p, a2, m])
        )
      )
    )
  );

  print("Primes in [5, 5000]: full search");
  printf("  Total (p,a₂) pairs with sf=-219: %d\n", count_small);
  printf("  Norm-form (p=k²+219) hits: %d\n", #norm_form_hits);
  printf("  Non-norm-form hits: %d\n", #non_norm_hits);
  print("");

  if(#norm_form_hits > 0,
    print("Norm-form hits:");
    for(i=1,#norm_form_hits,
      h = norm_form_hits[i];
      printf("  p=%-6d a2=%-6d m=%-4d k=%d (p=k²+219 ✓)\n", h[1],h[2],h[3],h[4])
    )
  );
  print("");

  if(#non_norm_hits > 0,
    print("Non-norm-form hits (NEW!):");
    for(i=1,#non_norm_hits,
      h = non_norm_hits[i];
      D = h[2]^2 - 4*h[1]^2;
      m_val = sqrtint(D / target_sf);
      printf("  p=%-6d a2=%-6d m=%-4d D=%-10d (not norm-form!)\n",
             h[1], h[2], h[3], D)
    ),
    print("No non-norm-form hits in [5,5000] — consistent with norm-form exclusivity.")
  );
}

print("");
print("=== Extended search: p in [5000, 10000] (spot-check) ===");

{
  target_sf = -219;
  sqr219 = sqrt(219);
  count2 = 0;
  new_hits = List();

  forprime(p = 5000, 10000,
    mmax = floor(2*p / sqr219);
    for(m = 1, mmax,
      rhs = 4*p^2 - 219*m^2;
      if(rhs <= 0, break);
      a2 = sqrtint(rhs);
      if(a2^2 != rhs, next);
      D = a2^2 - 4*p^2;
      sf = core(D);
      if(sf == target_sf,
        count2++;
        listput(new_hits, [p, a2, m])
      )
    )
  );

  printf("  Total pairs with sf=-219 in [5000,10000]: %d\n", count2);
  if(#new_hits > 0,
    for(i=1,#new_hits, h=new_hits[i];
      printf("  p=%d a2=%d m=%d\n",h[1],h[2],h[3])),
    print("  None found.")
  );
}

print("");
print("=== Theoretical note ===");
print("sf(a₂²-4p²) = -219 requires: a₂² + 219m² = 4p² for some integer m.");
print("This is the norm-N equation in Q(sqrt(-219)): N(a₂/2 + m*sqrt(-219)/2) = p².");
print("Solutions exist iff p splits or ramifies in Q(sqrt(-219)) = Q(sqrt(-3*73)).");
print("Since 219 = 3*73 and disc(Q(sqrt(-219))) = -219:");
print("  p splits in Q(sqrt(-219)) iff kronecker(-219, p) = 1.");
print("Count of solutions is related to h(-219) and the number of ideals of norm p.");
print("");

{
  \\ Compute h(-219) via bnfinit
  K = bnfinit(x^2 + 219);
  h = K.no;
  printf("Class number h(-219) = %d\n", h);
  printf("If h(-219) = 1, every split prime p is representable: infinitely many (p,a₂) pairs exist!\n");
  if(h == 1,
    print("h(-219)=1: the search above should find pairs for ALL split primes p."),
    printf("h(-219)=%d: representation is a coset condition; density 1/h among split primes.\n", h)
  );
}

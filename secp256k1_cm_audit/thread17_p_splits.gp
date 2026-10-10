\\ thread17_p_splits.gp
\\ Thread 17: Verify p always SPLITS in Q(sqrt(sf)) when p prime, a2≠0, p∤a2.
\\ Run: gp -q thread17_p_splits.gp

default(parisize, 64000000);

sf_part(n) = sign(n)*core(abs(n));

check_kron(p, a2) = {
  my(D = a2^2 - 4*p^2, sf, k);
  if(D == 0 || a2 % p == 0, return(-99));
  sf = sf_part(D);
  kronecker(sf, p)
};

print("=== Thread 17: p always SPLITS in Q(sqrt(sf)) when p∤a2 ===");

\\ --- Part A: 10 named cases ---
print("--- Part A: 10 named cases ---");

{
  my(tab, pass, fail, k);
  \\ [p, a2]
  tab = [[7,3],[11,4],[23,8],[47,22],[53,10],[97,12],[101,14],[199,8],[307,35],[2017,100]];
  pass = 0; fail = 0;
  for(i=1, #tab,
    k = check_kron(tab[i][1], tab[i][2]);
    printf("  p=%5d a2=%4d  kron=%2d  => %s\n",
      tab[i][1], tab[i][2], k, if(k==1,"SPLIT",if(k==0,"RAMIFIED","INERT or skip")));
    if(k==1, pass++, fail++);
  );
  printf("Part A: %d SPLIT, %d failures\n\n", pass, fail);
}

\\ --- Part B: Mass sweep 50 primes x up to 8 a2 values ---
print("--- Part B: Mass sweep ---");

{
  my(pass, fail, tot, p, a2, k, cnt);
  pass = 0; fail = 0; tot = 0; cnt = 0; p = 7;
  while(cnt < 50,
    if(isprime(p),
      cnt++;
      forstep(a2=1, min(p-1, 8), 1,
        if(a2 % p != 0,
          k = check_kron(p, a2);
          if(k != -99,
            tot++;
            if(k==1, pass++,
              printf("  FAIL: p=%d a2=%d kron=%d\n", p, a2, k);
              fail++;
            );
          );
        );
      );
    );
    p = nextprime(p+1);
  );
  printf("Part B: %d SPLIT, %d fail, %d total\n\n", pass, fail, tot);
}

\\ --- Part C: Ramification algebraic check ---
\\ p|sf requires p|D=a2^2-4p^2. Since p|4p^2, need p|a2^2 => p|a2. QED.
print("--- Part C: ramification algebraic check (verify p ∤ sf for all a2 in [1,p-1]) ---");

{
  my(found, sf);
  found = 0;
  forprime(p=5, 300,
    forstep(a2=1, p-1, 1,
      sf = sf_part(a2^2 - 4*p^2);
      if(sf % p == 0,
        printf("  UNEXPECTED ramification: p=%d a2=%d\n", p, a2);
        found = 1;
      );
    );
  );
  if(!found,
    print("  OK: no p|sf case found for p in [5,300], a2 in [1,p-1]");
  );
  print("");
}

\\ --- Part D: Inertness check (kron(sf,p)=-1?) ---
print("--- Part D: inertness check (all pairs from Part B should have kron=+1) ---");

{
  my(found, k, p, a2, cnt);
  found = 0; cnt = 0; p = 5;
  while(cnt < 100,
    if(isprime(p),
      cnt++;
      forstep(a2=1, min(p-1, 20), 1,
        if(a2 % p != 0,
          k = check_kron(p, a2);
          if(k == -1,
            printf("  UNEXPECTED inert: p=%d a2=%d\n", p, a2);
            found = 1;
          );
        );
      );
    );
    p = nextprime(p+1);
  );
  if(!found,
    print("  OK: no inertness found for p in [5,541], a2 in [1,min(p-1,20)]");
  );
}

print("\n=== CONCLUSION ===");
print("Algebraically: p|sf=>p|a2 (ramification impossible); inertness impossible by");
print("  beta-norm argument. Therefore p ALWAYS SPLITS in Q(sqrt(sf)) when p∤a2.");
print("Numerically verified above.");

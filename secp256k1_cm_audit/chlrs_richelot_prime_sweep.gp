\\ chlrs_richelot_prime_sweep.gp
\\
\\ Follow-up to chlrs_forward_map.gp's 0/12 root-choice falsification.
\\ Question: was howe_richelot_v5.gp / Test1's p=43 success (Richelot(alpha,beta)
\\ with alpha=cube root of -b1, beta=d*alpha, d=first QNR) a coincidence of small
\\ numbers, or does the naive construction actually work for a real subset of
\\ primes? Sweep small primes p=1 mod 3, several b1, and record the success rate.

default(parisize, 256000000);
read("howe_5pairs_v2.gp");

{
ntrials = 0; nsucc = 0; succlist = [];
forprime(pr = 11, 3000,
  if(pr%3 == 1,
    z3 = lift(polrootsmod(x^2+x+1, pr)[1]);
    for(b1 = 1, 7,
      if(b1 % pr != 0,
        dd = 2; while(kronecker(dd, pr) != -1 && dd < pr, dd++);
        if(dd < pr,
          b2 = lift(Mod(dd, pr)^3 * Mod(b1, pr));
          if(b2 != 0,
            rr = (-b1) % pr;
            sv = [0, (1+dd)%pr, 0]; qv = [0, 0, dd%pr];
            res = richelot(sv, qv, z3, rr, pr);
            if(res[1] != -1,
              E1 = ellinit([0, b1], pr);
              t1 = pr+1-ellcard(E1);
              target = (pr+1-t1)*(pr+1+t1);
              hh = Mod(1,pr)*x^6 + Mod(res[1],pr)*x^3 + Mod(res[2],pr);
              cp = hyperellcharpoly(hh);
              nj = subst(cp, variable(cp), 1);
              ntrials++;
              if(nj == target,
                nsucc++;
                succlist = concat(succlist, [[pr, b1, dd]]);
              );
            )
          )
        )
      )
    )
  )
);
print("Swept primes p=1 mod 3 in [11,3000], b1 in [1,7], d=first QNR mod p.");
print("trials=", ntrials, "  successes=", nsucc);
print("successes (p,b1,d): ", succlist);
}
quit;

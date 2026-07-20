\\ Thread 17: "p always splits in Q(sqrt(sf)) when p ndivides a2"
\\ Verifies 10 fresh (p,a2) + 3 large-p cases.
\\ Run: gp -q thread17_splits_verify.gp

default(parisize, 256000000);
default(timer, 0);

print("============================================================");
print("Thread 17: p-splits verification");
print("============================================================");
print("");

\\ squarefree_part(D) -> [sf, m] with D = sf * m^2, sf squarefree
squarefree_part(D) = {
  local(sf, m, D2, f, p2, e2);
  if (D == 0, error("squarefree_part: D=0"));
  D2 = abs(D);
  f  = factor(D2);
  sf = D2;
  m  = 1;
  forstep (k = 1, #f[,1], 1,
    p2 = f[k, 1];
    e2 = f[k, 2] \ 2;
    if (e2 > 0,
      m  = m * p2^e2;
      sf = sf / p2^(2*e2)
    )
  );
  if (D < 0, sf = -sf);
  return([sf, m]);
}

\\ disc_K(sf) = discriminant of Q(sqrt(sf))
disc_K(sf) = if (sf % 4 == 1, sf, 4*sf);

\\ check_one(p0, a0): returns [sf, kron, [P]^2_is_principal, ok]
check_one(p0, a0) = {
  local(D, sfm, sf, m, disc, kr, bnf0, fac_p, P, P2, res2, isprinc, ok);
  D    = a0^2 - 4*p0^2;
  sfm  = squarefree_part(D);
  sf   = sfm[1];
  m    = sfm[2];
  disc = disc_K(sf);
  kr   = kronecker(disc, p0);   \\ +1=splits, -1=inert, 0=ramified
  bnf0 = bnfinit(x^2 - sf, 1);
  fac_p = idealprimedec(bnf0, p0);
  P    = fac_p[1];
  P2   = idealpow(bnf0, P, 2);
  res2 = bnfisprincipal(bnf0, P2);
  isprinc = (norml2(res2[1]) == 0);
  ok = (kr == 1) && isprinc;
  return([sf, kr, isprinc, ok]);
}

\\ ---- 10 fresh test cases (p >= 307) ----
print("---- 10 fresh test cases ----");
print("  k | p   | a2   | sf    | kron | princ | ok");
print("  --|-----|------|-------|------|-------|----");

pass10 = 0; fail10 = 0;

run_case(idx, p0, a0) = {
  local(res, sf0, kr0, pr0, ok0, tag);
  if (a0 % p0 == 0, print("  ", idx, " | SKIP p|a2"); return(0));
  res = check_one(p0, a0);
  sf0 = res[1]; kr0 = res[2]; pr0 = res[3]; ok0 = res[4];
  tag = if (ok0, "PASS", "FAIL");
  print("  ", idx, " | ", p0, " | ", a0, " | ", sf0,
        " | ", kr0, " | ", if(pr0, 1, 0), " | ", tag);
  return(ok0);
}

pass10 += run_case(1,  307,   60);
pass10 += run_case(2,  311,   88);
pass10 += run_case(3,  313,  -44);
pass10 += run_case(4,  317,  100);
pass10 += run_case(5,  331,  -90);
pass10 += run_case(6,  337,   72);
pass10 += run_case(7,  347,  120);
pass10 += run_case(8,  349,  -66);
pass10 += run_case(9,  353,   40);
pass10 += run_case(10, 359, -160);

fail10 = 10 - pass10;
print("  => ", pass10, "/10 passed, ", fail10, " failed");
print("");

\\ ---- 3 large-p spot checks (p ~ 10^5) ----
print("---- Large-p spot check ----");

passL = 0;

run_large(p0, a0) = {
  local(res, ok0, tag);
  if (!isprime(p0), print("SKIP not prime: ", p0); return(0));
  if (a0 % p0 == 0, print("SKIP p|a2 for p=", p0); return(0));
  res = check_one(p0, a0);
  ok0 = res[4];
  tag = if (ok0, "PASS", "FAIL");
  print("  p=", p0, " a2=", a0, " sf=", res[1], " kron=", res[2],
        " princ=", if(res[3], 1, 0), " => ", tag);
  return(ok0);
}

passL += run_large(100003,   999);
passL += run_large(100019, -12345);
passL += run_large(100043,  7777);

failL = 3 - passL;
print("  => ", passL, "/3 passed, ", failL, " failed");
print("");

\\ ---- Algebraic summary ----
print("---- Algebraic argument (why p must split) ----");
print("(RAM) p|disc(K) iff p|sf; sf|D=a2^2-4p^2=a2^2 (mod p),");
print("      so p|sf => p|a2, excluded. Hence p NOT ramified.");
print("(INE) If p inert, only norm-p^2 ideal is (p); but (beta)!=(p)");
print("      since beta/p in O_K => p|a2, excluded. Hence NOT inert.");
print("=> p SPLITS. [P]^2=1 in Cl(K). QED.");
print("");

\\ ---- Final tally ----
total_new  = 13;
total_pass = pass10 + passL;
total_fail = fail10 + failL;

print("============================================================");
print("RESULT: ", total_pass, "/", total_new, " new cases passed.");
print("Combined with 82 prior cases: ", 82 + total_pass, " total.");
if (total_fail == 0, print("Zero failures. Theorem holds for all tested cases."));
if (total_fail > 0, print("ALERT: failures=", total_fail, ". Investigate."));
print("============================================================");
print("DONE.");

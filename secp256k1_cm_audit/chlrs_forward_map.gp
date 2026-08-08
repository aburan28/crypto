\\ chlrs_forward_map.gp
\\
\\ Thread 3 (CHLRS Igusa forward map) — proposed by RESEARCH_AUTOLAB_LOG.md 2026-07-27.
\\ Goal: given E1: y^2=x^3+b1, E2: y^2=x^3+b2 (j=0, non-degenerate pair, i.e.
\\ b2/b1 not a 6th power), find the correct branch-point identification (alpha,beta)
\\ so that the Z/3Z Richelot construction (howe_5pairs_v2.gp's richelot()) returns
\\ a genus-2 cover whose Jacobian matches E1 x E2.
\\
\\ Prior state (2026-07-27 log): howe_5pairs_v2.gp Test 2 (p=1009, d=11) used
\\ alpha^3=-b1, beta=d*alpha (so beta^3=-b2) and got a WRONG match (#Jac=1106283
\\ vs target 1018251). The log concluded a full CHLRS/Igusa-inversion formula was
\\ needed for the (r1,r2) branch parameters.
\\
\\ This script tests a cheaper hypothesis first: that the (alpha,beta) ROOT CHOICE
\\ among the 3 cube roots of (-b1) and (-b2) matters, i.e. richelot(sv,qv) is only
\\ correct for a specific relative "twist" between alpha and beta, not the naive
\\ one used in Test 2. We brute-force all 3 relative choices (beta = d*z3^k*alpha,
\\ k=0,1,2) and, for completeness, both sign conventions on rr (alpha^3 = +-b1).
\\ If none of these 6 combinations match, the root-choice hypothesis is falsified
\\ and the CHLRS branch-parameter formula (full paper apparatus) really is required.
\\
\\ Run: gp -q chlrs_forward_map.gp

default(parisize, 256000000);
default(timer, 0);
read("howe_5pairs_v2.gp"); \\ pulls in f3*, richelot(), check_jac() — also re-runs its own tests, harmless

print("");
print("================================================================");
print("chlrs_forward_map.gp: root-choice sweep for p=1009, pair (b1=11, d=11)");
print("================================================================");
{
  pp = 1009; b1 = 11; d = 11;      \\ same pair as Test 2 in howe_5pairs_v2.gp
  b2 = lift(Mod(d^3*b1, pp));
  z3roots = polrootsmod(x^2+x+1, pp);
  E1 = ellinit([0,b1], pp); E2 = ellinit([0,b2], pp);
  t1 = pp+1-ellcard(E1); t2 = pp+1-ellcard(E2);
  print("  b1=", b1, "  b2=", b2, "  trace(E1)=", t1, "  trace(E2)=", t2);
  target = (pp+1-t1)*(pp+1-t2);
  print("  target #Jac = (p+1-t1)(p+1-t2) = ", target);
  print("  z3 candidates: ", [lift(z3roots[1]), lift(z3roots[2])]);
  print("");

  nfound = 0;
  for(zi=1, 2,
    z3 = lift(z3roots[zi]);
    for(eps=0, 1,
      rr = if(eps==0, (-b1)%pp, b1%pp);
      for(k=0, 2,
        dk = lift(Mod(d,pp)*Mod(z3,pp)^k);
        sv = [0, (1+dk)%pp, 0];
        qv = [0, 0, dk%pp];
        res = richelot(sv, qv, z3, rr, pp);
        status = "n/a (cover not over F_p)";
        matchflag = 0;
        if(res[1]!=-1,
          chk = check_jac(res[1], res[2], t1, pp); \\ only chk[2] (#Jac) is used; chk[3] target is for the twist case, ignored
          nj = chk[2];
          matchflag = (nj == target);
          status = Str("aa=",res[1]," bb=",res[2]," #Jac=",nj," match=",matchflag);
        );
        print("  z3=",z3," eps=",eps," (rr=",if(eps==0,"-b1","+b1"),")  k=",k,"  beta=d*z3^",k,"*alpha:  ",status);
        if(matchflag, nfound++);
      )
    )
  );
  print("");
  if(nfound>0,
    print("  RESULT: ", nfound, "/12 root-choice combinations MATCH. Root-choice hypothesis CONFIRMED.");
  ,
    print("  RESULT: 0/12 root-choice combinations match target #Jac=", target, ".");
    print("  Root-choice hypothesis FALSIFIED. The relative alpha/beta twist is not");
    print("  the missing ingredient; a genuine CHLRS branch-parameter (r1,r2) formula,");
    print("  distinct from the naive cube roots of (-b1,-b2), is required.");
  );
}
print("");
print("================================================================");
print("Done.");
print("================================================================");

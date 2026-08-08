\\ chlrs_naive_cover_check.gp
\\
\\ Sanity check upstream of chlrs_forward_map.gp: does the NAIVE cover
\\ C0: y^2=(x^3+b1)(x^3+b2) itself (zero Richelot steps) have Jacobian
\\ L-polynomial equal to that of E1 x E2, for the same toy pairs used in
\\ howe_5pairs_v2.gp (p=43,b1=7,b2=13 and p=1009,b1=11,b2=515)?
\\
\\ Motivation: howe_richelot_v5.gp computes sigma_0/1/2, the 3 Richelot-DUAL
\\ neighbors of C0 in the (2,2)-isogeny graph. Since Richelot isogenies preserve
\\ the L-polynomial, sigma_0/1/2 must all share C0's own L-polynomial. If C0
\\ itself doesn't match target, none of its immediate Richelot neighbors can
\\ either -- which is exactly the 0/12 result chlrs_forward_map.gp found. This
\\ script checks the more basic question directly: does C0 match at all?

default(parisize, 256000000);

check_pair(pp, b1, b2, label) = {
  my(E1,E2,t1,t2,target,hh,cp,nj);
  E1 = ellinit([0,b1], pp); E2 = ellinit([0,b2], pp);
  t1 = pp+1-ellcard(E1); t2 = pp+1-ellcard(E2);
  target = (pp+1-t1)*(pp+1-t2);
  hh = Mod(1,pp)*x^6 + Mod(0,pp)*x^5 + Mod((b1+b2)%pp,pp)*x^3 + Mod((b1*b2)%pp,pp);
  \\ (x^3+b1)(x^3+b2) = x^6 + (b1+b2)x^3 + b1*b2
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  print("  [",label,"] p=",pp," b1=",b1," b2=",b2,"  t1=",t1," t2=",t2);
  print("     naive cover C0: y^2=x^6+",(b1+b2)%pp,"*x^3+",(b1*b2)%pp);
  print("     charpoly(C0) = ", cp);
  print("     #Jac(C0) = ", nj, "   target #E1*#E2 = ", target,
        "   match = ", nj==target);
  print("     factor over Q: ", factor(cp));
  print("");
};

print("================================================================");
print("Does the naive cover y^2=(x^3+b1)(x^3+b2) itself match Jac ~ E1 x E2?");
print("================================================================");
check_pair(43, 7, 13, "p=43 reference pair");
check_pair(1009, 11, 515, "p=1009 QNR-twist pair");

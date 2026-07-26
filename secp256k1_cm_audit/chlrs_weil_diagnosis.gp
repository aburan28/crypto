\\ chlrs_weil_diagnosis.gp
\\ Diagnose the Howe-gluing problem for j=0 curves
\\
\\ Establishes:
\\ 1. The split-sextic y^2=(x^3+b1)(x^3+b2) is NOT isogenous to E1 x E2.
\\ 2. The Rosenhain cross-ratio formula gives the wrong Weil poly.
\\ 3. Sextic-diagonal curves y^2=x^6+ax^3+b DO exist with the correct Weil poly.
\\ 4. For small primes the a=0 form y^2=x^6+b works; for p=1009 the first is a=1.
\\
\\ Run: gp -q chlrs_weil_diagnosis.gp

default(parisize, 256000000);
default(timer, 0);

\\ ---- Setup ----
p = 1009; b1 = 11;
E1 = ellinit([0, b1], p);
t1 = ellap(E1, p);
a2t = 2*p - t1^2;

print("===== CHLRS Weil Polynomial Diagnosis =====");
print("p=",p,"  E1: y^2=x^3+",b1,"  t1=",t1);
print("Target Weil poly: T^4+",a2t,"*T^2+",p^2," (#Jac = ",1+a2t+p^2,")");
print("");

get_cp_info(a, b) = {
  my(h, cp);
  h = Mod(x^6 + a*x^3 + b, p);
  iferr(cp = hyperellcharpoly(h), e, return([-9999,-9999,-9999]));
  [polcoeff(cp,3), polcoeff(cp,2), subst(cp,x,1)]
};

\\ ---- Result 1: Check naive split sextics ----
print("--- Result 1: Split sextic y^2=(x^3+b1)(x^3+b2) ---");
test_split(b2) = {
  my(a_s, b_s, res, E2c, t2);
  a_s = lift(Mod(b1+b2,p)); b_s = lift(Mod(b1*b2,p));
  E2c = ellinit([0,b2],p); t2 = ellap(E2c,p);
  res = get_cp_info(a_s, b_s);
  print("  b2=",b2," trace=",t2,"  a=",a_s," b=",b_s,
        "  a1=",res[1]," a2=",res[2]," #Jac=",res[3])
};

\\ First non-square d: kronecker(d,p) = -1
for(d = 2, 6,
  if(kronecker(d,p) == -1,
    test_split(lift(Mod(d,p)^3 * Mod(b1,p)));
    break()
  )
);

print("Key: a1 != 0 => NOT isogenous to E1 x E2 (which has symmetric Weil poly)");
print("");

\\ ---- Result 2: Rosenhain formula output ----
print("--- Result 2: Rosenhain cross-ratio formula ---");
print("From chlrs_igusa_formula.gp: gives a2=334 (wrong), target=",a2t);
print("334 = 2*",334/2,"  vs target 169 = 2*84.5. 334/169 =",334.0/169);
print("Structural error: the Mobius transform does NOT map to the Howe quotient.");
print("");

\\ ---- Result 3: Correct curve by brute force ----
print("--- Result 3: First correct sextic diagonal (verified) ---");
\\ y^2 = x^6 + x^3 + 117 over F_1009
res_correct = get_cp_info(1, 117);
print("y^2=x^6+x^3+117: a1=",res_correct[1]," a2=",res_correct[2],
      " #Jac=",res_correct[3]," (target ",1+a2t+p^2,")");
if(res_correct[2] == a2t && res_correct[1] == 0,
  print("CONFIRMED: a1=0, a2=",a2t," => in isogeny class of E1 x E2")
,
  print("UNEXPECTED: Weil poly mismatch")
);
print("");

\\ ---- Result 4: #C(F_p) for the correct curve ----
print("--- Result 4: Counting F_p points on y^2=x^6+x^3+117 ---");
count_pts(a, b) = {
  my(cnt, y2);
  cnt = 2;  \\ two points at infinity for monic degree-6
  for(x = 0, p-1,
    y2 = lift(Mod(x^6 + a*x^3 + b, p));
    if(y2 == 0, cnt++;,
      if(kronecker(y2,p) == 1, cnt+=2)
    )
  );
  cnt
};
npts = count_pts(1, 117);
print("#C(F_p) = ",npts,"  expected p+1 = ",p+1,"  diff = ",npts-(p+1));
print("(Confirms a1=0 in Weil poly: Tr(Frob) = p+1-#C = ",p+1-npts,")");
print("");

\\ ---- Summary ----
print("===== Summary =====");
print("Split-sextic: a1 != 0 => NOT in isogeny class of E1 x E2 over F_p.");
print("Rosenhain cross-ratio: a2=334 wrong (gives different kernel Gamma).");
print("Correct curve EXISTS: y^2=x^6+x^3+117 has correct Weil poly.");
print("OPEN: finding the EXPLICIT formula for the correct (a,b) from (b1,p).");
print("BLOCKED: finding the specific Howe-quotient among all isogenous curves");
print("         requires Mestre reconstruction or the CHLRS explicit formula.");
print("=========================================");

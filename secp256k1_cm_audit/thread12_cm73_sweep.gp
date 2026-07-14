\\ thread12_cm73_sweep.gp -- Thread 12: CM-73 extended sweep + CM field analysis
\\ Usage: gp --stacksize 32000000 -q thread12_cm73_sweep.gp

n73_count=0; nnf_count=0;

\\ Process each prime: print results inline, tally counts.
proc_p(p)={
  my(g,b1,b2,fld,P,f,a2,delta,rem,k2,k);
  if(p%6!=1,return(0));
  g=znprimroot(p); b1=lift(g^1); b2=lift(g^2);
  fld=ffgen(p,'t);
  P=fld^0*x^6+(b1+b2)*fld^0*x^3+(b1*b2%p)*fld^0;
  f=hyperellcharpoly(P);
  a2=polcoeff(f,2);
  delta=a2-2*p;
  rem=4*p-73;
  if(rem>0&&rem%3==0, k2=rem/3; k=sqrtint(k2); if(k^2==k2, nnf_count++; print("NF  p=",p,"  k=",k,"  a2-2p=",delta)));
  if(delta==-73, n73_count++; print("CM73  p=",p,"  4p-73=",4*p-73,"  kron(-219,p)=",kronecker(-219,p)));
  return(1);
}

print("=== Thread 12: CM-73 Extended Sweep p<=5000 ===");
print();
print("--- Part A+B: Sweep results (CM73=CM-73 hit, NF=norm-form hit) ---");
forprime(p=7,5000,proc_p(p));
print();
print("  CM-73 primes found: ",n73_count);
print("  Norm-form (4p=73+3k^2) primes found: ",nnf_count);
print();

\\ --- Part C: Class groups ---
print("--- Part C: Class groups ---");
K219=bnfinit(x^2+219,1); print("  Q(sqrt(-219)): h=",K219.no,"  clgp=",K219.clgp);
K73=bnfinit(x^2+73,1);   print("  Q(sqrt(-73)):  h=",K73.no,"  clgp=",K73.clgp);
K3=bnfinit(x^2+3,1);     print("  Q(sqrt(-3)):   h=",K3.no,"  clgp=",K3.clgp);
print();

\\ --- Part D: Kronecker for the 4 known CM-73 primes ---
print("--- Part D: Splitting analysis for CM-73 primes {19,37,79,109} ---");
print("  p=19:  kron(-219,19)=",kronecker(-219,19),"  kron(-73,19)=",kronecker(-73,19),"  kron(-3,19)=",kronecker(-3,19));
print("  p=37:  kron(-219,37)=",kronecker(-219,37),"  kron(-73,37)=",kronecker(-73,37),"  kron(-3,37)=",kronecker(-3,37));
print("  p=79:  kron(-219,79)=",kronecker(-219,79),"  kron(-73,79)=",kronecker(-73,79),"  kron(-3,79)=",kronecker(-3,79));
print("  p=109: kron(-219,109)=",kronecker(-219,109),"  kron(-73,109)=",kronecker(-73,109),"  kron(-3,109)=",kronecker(-3,109));
print("  -- Some norm-form primes (not CM-73) for comparison:");
print("  p=349: kron(-219,349)=",kronecker(-219,349),"  kron(-73,349)=",kronecker(-73,349));
print("  p=601: kron(-219,601)=",kronecker(-219,601),"  kron(-73,601)=",kronecker(-73,601));
print("  p=907: kron(-219,907)=",kronecker(-219,907),"  kron(-73,907)=",kronecker(-73,907));
print();

\\ --- Part E: Galois groups for CM-73 Weil polynomials ---
print("--- Part E: Galois groups for CM-73 Weil polys T^4+(2p-73)T^2+p^2 ---");
print("  p=19:  Gal=",polgalois(x^4+(2*19-73)*x^2+19^2));
print("  p=37:  Gal=",polgalois(x^4+(2*37-73)*x^2+37^2));
print("  p=79:  Gal=",polgalois(x^4+(2*79-73)*x^2+79^2));
print("  p=109: Gal=",polgalois(x^4+(2*109-73)*x^2+109^2));
print();

\\ --- Part F: Factor Weil poly over Q(sqrt(-219)) ---
print("--- Part F: Factoring Weil polys over Q(sqrt(-219)) ---");
\\ Use 'y' for the field generator to avoid variable collision with 'x'
K=nfinit(y^2+219);
print("  p=19:  T^4-35T^2+361  factors=",nffactor(K,x^4-35*x^2+361));
print("  p=37:  T^4+T^2+1369   factors=",nffactor(K,x^4+x^2+1369));
print("  p=79:  T^4+85T^2+6241 factors=",nffactor(K,x^4+85*x^2+6241));
print("  p=109: T^4+145T^2+11881 factors=",nffactor(K,x^4+145*x^2+11881));
print();

\\ --- Part G: Subfield structure of the quartic CM splitting field ---
\\ The splitting field L of T^4+(2p-73)T^2+p^2 over Q has Gal(L/Q)=V_4.
\\ Its three quadratic subfields for p=19 (k=1):
\\   F1 = Q(sqrt(-219))  -- from alpha - alpha_bar = sqrt(-219)
\\   F2 = Q(sqrt(73))    -- from (beta + p/beta)^2 = 73  where beta=sqrt(alpha)
\\   F3 = Q(sqrt(-3))    -- from (beta - p/beta)^2 = -3
\\ Verify: (-219) = (-3)*73, so F1 = F2.F3 (compositum). Confirmed.
print("--- Part G: Subfield structure verification ---");
print("  (beta + p/beta)^2 for p=19,k=1: alpha=(35+sqrt(-219))/2");
print("  = alpha + 2p + p^2/alpha = alpha + 2p + alpha_bar = 2*Re(alpha) + 2p");
print("  = 2*(35/2) + 2*19 = 35 + 38 = 73. So F2 = Q(sqrt(73)). Check.");
print("  (beta - p/beta)^2 = alpha - 2p + alpha_bar = (alpha+alpha_bar) - 2p = 35 - 38 = -3.");
print("  So F3 = Q(sqrt(-3)). Check.");
print("  Compositum F1 = F2*F3 = Q(sqrt(73), sqrt(-3)) has product sqrt(-219). Check.");
print();

print("=== Thread 12 complete ===");
quit();

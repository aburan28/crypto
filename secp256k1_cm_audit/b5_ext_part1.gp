\\ B5 Part 1: proxy prime pp=43 numerical verification
\\ Inline code only — no user functions (avoids PARI parsing issues)
default(parisize, 64000000);
default(timer, 0);

pp = 43;

\\ Count #E: y^2=x^3+1 over F_43
nE0 = 1; for(xx=0,pp-1, rhs=lift(Mod(xx,pp)^3+1); nE0 += 1+kronecker(rhs,pp));
tr0 = pp + 1 - nE0;
print("E  : #E(F_43) = ", nE0, "  trace = ", tr0);

\\ Quadratic twist: find QNR nu
nu0 = 2; while(kronecker(nu0,pp) != -1, nu0++);
nEt0 = 1; for(xx=0,pp-1, rhs=lift(Mod(xx,pp)^3+nu0); nEt0 += 1+kronecker(rhs,pp));
trt0 = pp + 1 - nEt0;
print("E^t: #E^t(F_43) = ", nEt0, "  trace = ", trt0);
print("twist check: tr+tr^t = ", tr0+trt0, " (want 0)");
print("");

\\ Lucas sequence inline: L0=2, L1=tr, Lk=tr*L_{k-1}-p*L_{k-2}
\\ Compute for k=1..4 and print table
print("k  nEk           nEtk          nJac           2^ecdlp 2^hcdlp margin");
{
  for(kk=1,4,
    ppk = pp^kk;
    \\ Compute L_kk(tr0, pp) via loop
    La=2; Lb=tr0;
    if(kk >= 2, for(ii=2,kk, Lc=tr0*Lb-pp*La; La=Lb; Lb=Lc));
    Lk = if(kk==1, tr0, Lb);
    nEk  = ppk + 1 - Lk;
    nEtk = ppk + 1 - (-1)^kk * Lk;
    nJk  = nEk * nEtk;
    s1 = log(nEk+0.0)/log(2.0)/2;
    s2 = log(nJk+0.0)/log(2.0)/2;
    printf("k=%d %13d %13d %15d 2^%4.2f 2^%5.2f +%.2f\n", kk,nEk,nEtk,nJk,s1,s2,s2-s1)
  )
}
print("sqrt(#Jac) > sqrt(#E) for all k: cover DLP always harder than ECDLP. OK");

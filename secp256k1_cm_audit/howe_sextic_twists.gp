\\ Thread 18: Howe (H1)+(H2)+(H3) pairwise check for ALL 15 pairs
\\ of the 6 sextic twists of secp256k1 (j=0, CM by Q(sqrt(-3))).
\\
\\ PARI/GP 2.15 NOTE: for/while bodies must be on ONE line (no multi-line).
\\ Run: gp -q howe_sextic_twists.gp

default(parisize, 128000000);
default(timer, 0);

p  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n0 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t0 = p + 1 - n0;
print("p = ", p); print("t_0 = ", t0);

\\ CM decomposition: 4p = t^2 + 3*b^2
disc4 = 4*p - t0^2;
b = sqrtint(disc4/3);
print("b = ", b, "  (CM param, t^2+3b^2==4p: ", t0^2+3*b^2==4*p, ")");

\\ 6 sextic-twist traces using individual assignments
traces = vector(6);
traces[1] = t0; traces[2] = (t0-3*b)/2; traces[3] = -(t0+3*b)/2;
traces[4] = -t0; traces[5] = (3*b-t0)/2; traces[6] = (t0+3*b)/2;
orders = vector(6, k, p+1-traces[k]);
print("orders[1]==n0: ", orders[1]==n0, "  orders[4]==p+1+t0: ", orders[4]==p+1+t0);
print("");
print("--- 6 sextic twist group orders ---");
for (k=1, 6, print("  E_",k-1,":  t=",traces[k],"  #E=",orders[k]));
print("");

\\ Find primitive 6th root of unity d mod p (single-line body required)
d6 = 0; kk = 2;
while (d6 == 0, dc = Mod(kk,p)^((p-1)/6); if (dc^3 != Mod(1,p) && dc^2 != Mod(1,p), d6 = lift(dc); print("d = ", kk, "^{(p-1)/6}  d^6==1: ", Mod(d6,p)^6==Mod(1,p))); kk = kk+1);
print("");

\\ 2-torsion cubic character: -B_k is a cube iff (-B_k)^{(p-1)/3} == 1 mod p
\\ Compute for each twist k=0..5 without a loop (to avoid multi-line body issues)
exp3 = (p-1)/3;
d6p = Mod(d6,p);
Bvec = vector(6); Bvec[1]=Mod(7,p); Bvec[2]=Mod(7,p)*d6p; Bvec[3]=Mod(7,p)*d6p^2; Bvec[4]=Mod(7,p)*d6p^3; Bvec[5]=Mod(7,p)*d6p^4; Bvec[6]=Mod(7,p)*d6p^5;
is_cube = vector(6);
for (k=1, 6, is_cube[k] = ((-Bvec[k])^exp3 == Mod(1,p)));
print("--- 2-torsion cubic character (-B_k is a cube?) ---");
for (k=1, 6, print("  E_",k-1,":  cube=",is_cube[k],"  x^3+B_",k-1," is ",if(is_cube[k],"SPLIT","IRREDUCIBLE")));
print("");

\\ All 15 pairwise Howe checks (nested for, single-line bodies)
print("=== All 15 pairwise Howe (H1)+(H2)+(H3) checks ===");
n_pass = 0; n_fail = 0;
for (i=0, 4, for (j=i+1, 5, ni=orders[i+1]; nj=orders[j+1]; H1=(ni!=nj); g12=gcd(ni,nj); H3=(g12==1); H2=(is_cube[i+1]==is_cube[j+1]); apass=H1&&H2&&H3; if(apass,n_pass=n_pass+1,n_fail=n_fail+1); print("(E_",i,",E_",j,"): H1=",H1," H2=",H2," H3(gcd=",g12,")=",H3,"  ALL=",apass)));
print("");
print("==> ", n_pass, " / 15 pairs satisfy all 3 Howe conditions; ", n_fail, " fail.");

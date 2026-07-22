\\ Thread 18 compact: compute all 15 pair Howe conditions.
\\ Uses {} blocks and single-line if() to avoid PARI parsing issues.

default(parisize, 512000000);

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
t0 = p + 1 - n;

\\ Eisenstein factorization
disc_eis = 4*p - t0^2;
u_eis = sqrtint(disc_eis\3);
a_eis = (t0 + u_eis)\2; b_eis = 2*a_eis - t0;
if(a_eis^2 - a_eis*b_eis + b_eis^2 != p, a_eis = (t0-u_eis)\2; b_eis = 2*a_eis - t0);

print("Eisenstein: a=", a_eis, " b=", b_eis);
print("Check: a^2-ab+b^2=p: ", a_eis^2-a_eis*b_eis+b_eis^2 == p, "  2a-b=t0: ", 2*a_eis-b_eis == t0);
print("");

\\ Six sextic twist traces: {+/-s1, +/-s2, +/-s3}
\\ with s1+s2+s3 related by: (2a-b) + (-(a+b)) + (2b-a) = 0
s1 = 2*a_eis - b_eis;
s2 = -(a_eis + b_eis);
s3 = 2*b_eis - a_eis;
print("s1=", s1);
print("s2=", s2);
print("s3=", s3);
print("s1+s2+s3=", s1+s2+s3, " (should be 0)");
print("");

\\ T[i] = trace of E_i, N[i] = #E_i = p+1-T[i]
\\ Labelling: 1=s1, 2=-s1, 3=s2, 4=-s2, 5=s3, 6=-s3
T = [s1, -s1, s2, -s2, s3, -s3];
N = vector(6, i, p+1-T[i]);

print("Six group orders:");
{
for(i=1,6,
  print("  E_",i,": #E=",N[i],"  trace=",T[i],"  4|trace=", T[i]%4==0));
}
print("");

\\ H2 proxy: 4|T[i] <=> split 2-torsion (since p=3 mod 4 so p+1=0 mod 4)
\\ Theory: x^3+c splits completely over F_p iff -c is a cube in F_p*.
\\ When x^3+c is irred, Galois acts as Z/3Z on 3 non-trivial 2-tors pts.
\\ When x^3+c splits, E[2](F_p) = Z/2Z x Z/2Z, so 4 | #E.
\\ Among 6 twists (p=1 mod 3 so F_p has cube roots of unity):
\\   exactly 2 have split 2-torsion (the pair +-s3 here since 4|s3),
\\   4 have irreducible 2-torsion.

\\ Check all 15 pairs
print("=================================================");
print("All 15 pairs (H1, H2, H3) check:");
print("H1: N[i]!=N[j]  H2: 4|T[i] iff 4|T[j]  H3: gcd=1");
print("(ss)=both split  (ii)=both irred  (si)=mixed");
print("=================================================");
print("");

{
for(i=1,5,
  for(j=i+1,6,
    h1 = (N[i] != N[j]);
    si_split = (T[i] % 4 == 0);
    sj_split = (T[j] % 4 == 0);
    h2 = (si_split == sj_split);
    g12 = gcd(N[i], N[j]);
    h3 = (g12 == 1);
    tag = "(si)";
    if(si_split && sj_split, tag = "(ss)");
    if(!si_split && !sj_split, tag = "(ii)");
    print("(",i,",",j,")",tag,"  H1=",if(h1==1,"Y","N")," H2=",if(h2==1,"Y","N")," H3=",if(h3==1,"Y","N"),"  gcd=",g12)));
}
print("");

\\ Summary counts
{
n_h2 = 0; n_all = 0;
for(i=1,5, for(j=i+1,6,
  si_split = (T[i]%4==0); sj_split = (T[j]%4==0); h2 = (si_split==sj_split);
  g12 = gcd(N[i],N[j]); h3 = (g12==1); h1 = (N[i]!=N[j]);
  if(h2, n_h2++);
  if(h1 && h2 && h3, n_all++)));
print("Pairs satisfying H2 only: ", n_h2, " / 15");
print("Pairs satisfying H1+H2+H3: ", n_all, " / 15");
}
print("");

\\ Factor the 4 irred-twist group orders (trial division to 10^6)
print("Factorisation (to 10^6 bound) of irred-twist group orders:");
{
for(i=1,6,
  if(T[i]%4 != 0,
    f = factorint(N[i], 10^6);
    print("  N[",i,"] = ", f)));
}
print("");

\\ Cross-check with howe_gluing_test.gp
print("Cross-check: gcd(n, n_twist) = gcd(N[1],N[2]) = ", gcd(N[1],N[2]));
print("(Should be 1 per howe_gluing_test.gp)");
print("Done.");

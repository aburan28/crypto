\\ howe_5pairs_v2.gp
\\
\\ FIXED version: pass prime 'pp' as explicit parameter to all F_{pp^3} functions.
\\ Investigates the 5 Howe-qualifying secp256k1 sextic-twist pairs.
\\
\\ Key results expected:
\\   - p=43 toy:   richelot gives a=41, b=5 (verified by howe_richelot_v5.gp)
\\   - p=1009 toy: richelot gives correct Howe cover for E1 x twist
\\   - Pair (0,3): degenerate (d=-1, sv=0, Delta=0); needs different approach
\\   - Pair (0,1)/(0,4): d is a cube root of h; need to check if it's a non-square
\\   - Naive cover: 189=27*7 — check which sextic twist class it falls in
\\
\\ Run: gp -q howe_5pairs_v2.gp

default(parisize, 256000000);
default(timer, 0);

\\ ================================================================
\\ F_{pp^3} arithmetic: prime pp, reduction a^3 = rr.
\\ Elements [c0,c1,c2] = c0 + c1*a + c2*a^2.
\\ ================================================================

f3add(u, v, pp) = [(u[1]+v[1])%pp, (u[2]+v[2])%pp, (u[3]+v[3])%pp];
f3neg(u, pp)    = [(-u[1])%pp, (-u[2])%pp, (-u[3])%pp];
f3scl(c, u, pp) = [(c*u[1])%pp, (c*u[2])%pp, (c*u[3])%pp];

f3mul(u, v, rr, pp) = {
  my(c0, c1, c2, c3, c4);
  c0 = (u[1]*v[1]) % pp;
  c1 = (u[1]*v[2] + u[2]*v[1]) % pp;
  c2 = (u[1]*v[3] + u[2]*v[2] + u[3]*v[1]) % pp;
  c3 = (u[2]*v[3] + u[3]*v[2]) % pp;
  c4 = (u[3]*v[3]) % pp;
  [(c0 + rr*c3) % pp,  (c1 + rr*c4) % pp,  c2 % pp]
};

f3inv(u, rr, pp) = {
  my(a, b, c, nrm, ni);
  a = u[1]; b = u[2]; c = u[3];
  \\ Norm of a+b*alpha+c*alpha^2 in F_pp[alpha]/(alpha^3-rr)
  \\ N = a^3 + rr*b^3 + rr^2*c^3 - 3*rr*a*b*c
  nrm = lift(Mod(a^3 + rr*b^3 + rr^2*c^3 - 3*rr*a*b*c, pp));
  if(nrm == 0, error("f3inv: norm=0, element not invertible"));
  ni = lift(Mod(nrm, pp)^(-1));
  [((a^2 - rr*b*c)*ni) % pp,
   ((rr*c^2 - a*b)*ni) % pp,
   ((b^2 - a*c)*ni) % pp]
};

\\ Z/3Z Richelot: given alpha+beta and alpha*beta in F_{pp^3},
\\ compute the Howe-glued dual cover y^2 = x^6 + aa*x^3 + bb.
\\ Returns [aa, bb] if the cover is over F_pp, else [-1, -1].

richelot(sv, qv, z3, rr, pp) = {
  my(z3sq, G1c, G1x, G2c, G2x, G3c, G3x, D0, D0inv);
  my(H1n2,H1n1,H1n0,H2n2,H2n1,H2n0,H3n2,H3n1,H3n0);
  my(H1x2,H1x1,H1x0,H2x2,H2x1,H2x0,H3x2,H3x1,H3x0);
  my(P0,P1,P2,P3,P4,Q0r,Q1r,Q2r,Q3r,Q4r,Q5r,Q6r,lc,lcinv,aa,bb);

  z3sq = lift(Mod(z3, pp)^2);

  G1c = qv;          G1x = f3neg(sv, pp);
  G2c = f3scl(z3sq, qv, pp);  G2x = f3neg(f3scl(z3, sv, pp), pp);
  G3c = f3scl(z3, qv, pp);   G3x = f3neg(f3scl(z3sq, sv, pp), pp);

  D0 = f3add(f3add(
    f3add(f3mul(G2x,G3c,rr,pp), f3neg(f3mul(G3x,G2c,rr,pp),pp), pp),
    f3neg(f3mul(G1x, f3add(G3c, f3neg(G2c,pp), pp), rr, pp), pp), pp),
    f3mul(G1c, f3add(G3x, f3neg(G2x,pp), pp), rr, pp), pp);

  D0inv = f3inv(D0, rr, pp);

  H1n2 = f3add(G3x, f3neg(G2x,pp), pp);
  H1n1 = f3scl(2, f3add(G3c, f3neg(G2c,pp), pp), pp);
  H1n0 = f3add(f3mul(G2x,G3c,rr,pp), f3neg(f3mul(G2c,G3x,rr,pp),pp), pp);
  H2n2 = f3add(G1x, f3neg(G3x,pp), pp);
  H2n1 = f3scl(2, f3add(G1c, f3neg(G3c,pp), pp), pp);
  H2n0 = f3add(f3mul(G3x,G1c,rr,pp), f3neg(f3mul(G3c,G1x,rr,pp),pp), pp);
  H3n2 = f3add(G2x, f3neg(G1x,pp), pp);
  H3n1 = f3scl(2, f3add(G2c, f3neg(G1c,pp), pp), pp);
  H3n0 = f3add(f3mul(G1x,G2c,rr,pp), f3neg(f3mul(G1c,G2x,rr,pp),pp), pp);

  H1x2=f3mul(H1n2,D0inv,rr,pp); H1x1=f3mul(H1n1,D0inv,rr,pp); H1x0=f3mul(H1n0,D0inv,rr,pp);
  H2x2=f3mul(H2n2,D0inv,rr,pp); H2x1=f3mul(H2n1,D0inv,rr,pp); H2x0=f3mul(H2n0,D0inv,rr,pp);
  H3x2=f3mul(H3n2,D0inv,rr,pp); H3x1=f3mul(H3n1,D0inv,rr,pp); H3x0=f3mul(H3n0,D0inv,rr,pp);

  P0=f3mul(H1x0,H2x0,rr,pp);
  P1=f3add(f3mul(H1x0,H2x1,rr,pp),f3mul(H1x1,H2x0,rr,pp),pp);
  P2=f3add(f3add(f3mul(H1x0,H2x2,rr,pp),f3mul(H1x1,H2x1,rr,pp),pp),f3mul(H1x2,H2x0,rr,pp),pp);
  P3=f3add(f3mul(H1x1,H2x2,rr,pp),f3mul(H1x2,H2x1,rr,pp),pp);
  P4=f3mul(H1x2,H2x2,rr,pp);

  Q0r=f3mul(P0,H3x0,rr,pp);
  Q1r=f3add(f3mul(P0,H3x1,rr,pp),f3mul(P1,H3x0,rr,pp),pp);
  Q2r=f3add(f3add(f3mul(P0,H3x2,rr,pp),f3mul(P1,H3x1,rr,pp),pp),f3mul(P2,H3x0,rr,pp),pp);
  Q3r=f3add(f3add(f3mul(P1,H3x2,rr,pp),f3mul(P2,H3x1,rr,pp),pp),f3mul(P3,H3x0,rr,pp),pp);
  Q4r=f3add(f3add(f3mul(P2,H3x2,rr,pp),f3mul(P3,H3x1,rr,pp),pp),f3mul(P4,H3x0,rr,pp),pp);
  Q5r=f3add(f3mul(P3,H3x2,rr,pp),f3mul(P4,H3x1,rr,pp),pp);
  Q6r=f3mul(P4,H3x2,rr,pp);

  if(Q0r[2]!=0||Q0r[3]!=0||Q3r[2]!=0||Q3r[3]!=0||Q6r[2]!=0||Q6r[3]!=0,return([-1,-1]));
  lc = Q6r[1]; if(lc==0, return([-1,-1]));
  lcinv = lift(Mod(lc,pp)^(-1));
  aa = (Q3r[1]*lcinv)%pp; bb=(Q0r[1]*lcinv)%pp;
  [aa, bb]
};

\\ Check: given (aa,bb), compute Frobenius poly and compare #Jac vs target.
check_jac(aa, bb, t_expected, pp) = {
  my(hh, cp, nj, target);
  hh = Mod(1,pp)*x^6 + Mod(aa,pp)*x^3 + Mod(bb,pp);
  cp = hyperellcharpoly(hh);
  nj = subst(cp, variable(cp), 1);
  target = (pp+1-t_expected) * (pp+1+t_expected);
  [cp, nj, target, nj==target]
};

\\ ================================================================
print("================================================================");
print("Test 1: p=43, E1: y^2=x^3+7, E2: y^2=x^3+13 (d=2)");
print("Expected: a=41, b=5  (from howe_richelot_v5.gp)");
print("================================================================");
{
  pp=43; rr=36; z3=6;  \\ a^3=36=-7 mod 43; cube root of unity z3=6
  sv=[0,3,0]; qv=[0,0,2];  \\ sv=(1+2)*alpha=3a, qv=2*alpha^2=2a^2
  res=richelot(sv,qv,z3,rr,pp);
  print("  Richelot(sv=[0,3,0],qv=[0,0,2]): a=",res[1],"  b=",res[2]);
  if(res[1]==41&&res[2]==5, print("  ✓ CORRECT"), print("  ✗ WRONG"));
  if(res[1]!=-1,
    chk=check_jac(res[1],res[2],13,pp);
    print("  Frobenius poly: ",chk[1]);
    print("  #Jac=",chk[2],"  target=",chk[3],"  match=",chk[4]);
  );
}
print("");

\\ ================================================================
print("================================================================");
print("Test 2: p=1009, E1: y^2=x^3+11, E2: y^2=x^3+515 (d=11)");
print("================================================================");
{
  pp=1009; b1=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  print("  z3 mod ",pp," = ",z3);
  d=2; while(kronecker(d,pp)!=-1, d++);
  print("  First non-square d = ",d);
  b2=lift(Mod(d^3*b1,pp));
  print("  E2 coefficient: ",b2);
  E1=ellinit([0,b1],pp); E2=ellinit([0,b2],pp);
  t1=pp+1-ellcard(E1); t2=pp+1-ellcard(E2);
  print("  trace E1=",t1,"  trace E2=",t2);
  rr=(-b1)%pp;
  sv=[0,(1+d)%pp,0]; qv=[0,0,d%pp];
  res=richelot(sv,qv,z3,rr,pp);
  print("  Richelot: a=",res[1],"  b=",res[2]);
  if(res[1]!=-1,
    chk=check_jac(res[1],res[2],t1,pp);
    print("  Frobenius poly: ",chk[1]);
    print("  #Jac=",chk[2],"  target=",chk[3],"  match=",chk[4]);
    if(chk[4], print("  ✓ CORRECT"), print("  ✗ WRONG"));
  );
}
print("");

\\ ================================================================
print("================================================================");
print("Test 3: Pair (0,3) degenerate case (d=-1, sv=0).");
print("  E0: y^2=x^3+7, E3: y^2=x^3-7 (quadratic twist via h^3=-1)");
print("  The standard Z/3Z Richelot degenerates; use alt Moebius form.");
print("================================================================");
print("");
print("  With d=-1: G_i = (x-alpha)(x+alpha) = x^2-alpha^2.");
print("  The Richelot image is NOT in the x^6+ax^3+b form.");
print("  For d=-1, use the factorization of (x^3+b)(x^3-b) = x^6-b^2.");
print("  This gives y^2 = x^6 - b^2 = (x^3-b)(x^3+b) — trivially split.");
print("  The Richelot of y^2=G1*G2*G3 with G_i=x^2-z3^{i-1}*alpha^2 is:");
print("  G1=(x-alpha)(x+alpha), G2=(x-z3*alpha)(x+z3*alpha), G3=(x-z3^2*alpha)(x+z3^2*alpha)");
print("  = x^2-alpha^2, x^2-z3^2*alpha^2, x^2-z3*alpha^2.");
print("  Product: (x^2-alpha^2)(x^2-z3^2*alpha^2)(x^2-z3*alpha^2) = x^6-alpha^6.");
print("  So y^2=x^6-b^2 has Richelot dual related to b^{1/3} factorization.");
print("");
print("  This is a DEGENERATE case: the standard Howe-gluing via (x^3+b)(x^3-b)");
print("  gives a SPLIT Jacobian (Jac(C) ~ E x E) only IF the curve is smooth");
print("  and the (2,2)-isogeny is non-trivial. But x^6-b^2 = (x^3+b)(x^3-b)");
print("  might have trivial Richelot (mapping to product of degree-1 factors).");
print("");
print("  For secp256k1 pair (0,3): b=7, h^3=-1, so the cover is y^2=x^6-49.");
print("  BLOCKED: This degenerate case needs the THETA FUNCTION or MESTRE approach");
print("  to find the correct non-trivial (2,2)-isogeny.");
print("  See RESEARCH_MESTRE_HOWE.md for alternative reconstruction algorithms.");
print("");

\\ ================================================================
print("================================================================");
print("Test 4: Pair (0,1) for p=1009");
print("  E0: y^2=x^3+11, E1: y^2=x^3+11*h for h=6th root of unity mod 1009");
print("================================================================");
{
  pp=1009; b0=11;
  z3=lift(polrootsmod(x^2+x+1,pp)[1]);
  h=lift(polrootsmod(x^2-x+1,pp)[1]);  \\ primitive 6th root
  print("  h = ",h,"  h^3 mod p = ",lift(Mod(h,pp)^3));
  b1=lift(Mod(b0*h,pp));
  print("  E1: b1 = 11*h mod 1009 = ",b1);

  \\ Check if h is a cube mod pp
  h_cube_test = lift(Mod(h,pp)^((pp-1)/3));
  print("  h^{(p-1)/3} mod p = ",h_cube_test,"  (=1 iff h is a cube)");

  if(h_cube_test==1,
    \\ Find cube root d of h
    for(d_cand=1, pp-1,
      if(lift(Mod(d_cand,pp)^3)==h%pp,
        print("  Cube root of h: d=",d_cand,"  kr(d,p)=",kronecker(d_cand,pp));
        \\ Try d, d*z3, d*z3^2 to find a non-square
        d0=d_cand;
        d1=lift(Mod(d_cand*z3,pp)); d2=lift(Mod(d_cand*z3^2,pp));
        print("  d*z3  = ",d1,"  kr=",kronecker(d1,pp));
        print("  d*z3^2= ",d2,"  kr=",kronecker(d2,pp));
        \\ Pick the non-square (if any)
        d_ns = 0;
        if(kronecker(d0,pp)==-1, d_ns=d0);
        if(kronecker(d1,pp)==-1, d_ns=d1);
        if(kronecker(d2,pp)==-1, d_ns=d2);
        if(d_ns>0,
          print("  Using d_ns=",d_ns,"  d_ns^3 mod p=",lift(Mod(d_ns,pp)^3),
                "  (expect h=",h%pp,")");
          rr01=(-b0)%pp;
          sv01=[0,(1+d_ns)%pp,0]; qv01=[0,0,d_ns%pp];
          res01=richelot(sv01,qv01,z3,rr01,pp);
          print("  Richelot(pair 0,1): a=",res01[1],"  b=",res01[2]);
          if(res01[1]!=-1,
            E0_tmp=ellinit([0,b0],pp); E1_tmp=ellinit([0,b1],pp);
            t0_tmp=pp+1-ellcard(E0_tmp); t1_tmp=pp+1-ellcard(E1_tmp);
            print("  t0=",t0_tmp,"  t1=",t1_tmp);
            chk01=check_jac(res01[1],res01[2],t0_tmp,pp);
            print("  Frobenius poly: ",chk01[1]);
            nj01=chk01[2]; target01=(pp+1-t0_tmp)*(pp+1-t1_tmp);
            print("  #Jac=",nj01,"  (p+1-t0)(p+1-t1)=",target01,"  match=",nj01==target01);
            \\ Also check (p+1-t0)(p+1+t0) in case it's isogenous to E0 x E0-twist
            print("  (p+1-t0)(p+1+t0)=",((pp+1-t0_tmp)*(pp+1+t0_tmp)));
            if(nj01==target01,print("  ✓ Jac ~ E0 x E1"),
              if(nj01==(pp+1-t0_tmp)*(pp+1+t0_tmp),print("  ? Jac ~ E0 x E0-twist"),
                print("  ? Jac ~ something else"))),
            print("  Cover not over F_p")
          ),
          print("  All cube roots of h are SQUARES mod p — pair (0,1) cover may not be over F_p")
        );
        break
      )
    ),
    print("  h is NOT a cube mod p; cube root in F_{p^3} only.")
  );
  print("");
}

\\ ================================================================
print("================================================================");
print("Test 5: Is 189 one of the 6 secp256k1 sextic twist classes?");
print("  Naive cover: y^2=(x^3+7)(x^3+189), d=3, b2=189=27*7.");
print("  Check: 189/c_k is a 6th power mod p for some k=0..5.");
print("================================================================");
{
  pp=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
  h=lift(polrootsmod(x^2-x+1,pp)[1]);
  b_naive=189;
  found=0;
  for(kk=0,5,
    c_k=lift(Mod(7,pp)*Mod(h,pp)^kk);
    ratio=lift(Mod(b_naive,pp)*Mod(c_k,pp)^(-1));
    test6=lift(Mod(ratio,pp)^((pp-1)/6));
    if(test6==1,
      print("  b_naive=189 is F_p-isomorphic to twist class k=",kk,
            " (since 189/c_k = 6th power, ratio=",ratio,")");
      found=1;
    );
  );
  if(!found,
    print("  189 is NOT F_p-isomorphic to any of the 6 canonical sextic twist classes.");
    print("  This means the naive cover y^2=(x^3+7)(x^3+189) uses a DIFFERENT quadratic");
    print("  twist from the 5 H-qualifying pairs in Thread 18.");
    print("  (Equivalently: 189/7=27=3^3 is not a 6th power times any h^k/h^0=h^k.)");
  );
  \\ Also check: what is the quadratic non-residue character of 3 and 27?
  print("");
  print("  kronecker(3, p_secp) = ",kronecker(3,pp));
  print("  kronecker(27, p_secp) = ",kronecker(27,pp)," = kronecker(3,p)^3");
  print("  (27 non-square iff 3 is non-square)");
  print("");
  \\ kronecker(3,p) = ?  p = 2^256-2^32-977
  \\ By quadratic reciprocity and p mod 12:
  print("  p mod 12 = ",pp%12);
  print("  p mod 3 = ",pp%3);
}
print("");

print("================================================================");
print("Done.");
print("================================================================");

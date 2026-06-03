\\ howe_richelot_v2.gp — flat single-line style for PARI 2.15.4
\\ Compute Richelot image of y^2=(x^3+7)(x^3+13) over F_{43^3}
\\ to identify the Howe-glued genus-2 curve.
\\ Run: gp -q howe_richelot_v2.gp

default(parisize, 64000000);
default(timer, 0);
p = 43;

\\ F_{43^3} = F_43[a]/(a^3 - 36). Elements: [c0,c1,c2] = c0+c1*a+c2*a^2.
\\ Helpers: defined as one-liners to avoid PARI multi-line parse bug.
fadd(u,v) = [(u[1]+v[1])%p, (u[2]+v[2])%p, (u[3]+v[3])%p];
fneg(u) = [(-u[1])%p, (-u[2])%p, (-u[3])%p];
fscl(c,u) = [(c*u[1])%p, (c*u[2])%p, (c*u[3])%p];
fcon(c) = [c%p, 0, 0];
fmul(u,v) = [(u[1]*v[1]+36*(u[2]*v[3]+u[3]*v[2]))%p, (u[1]*v[2]+u[2]*v[1]+36*u[3]*v[3])%p, (u[1]*v[3]+u[2]*v[2]+u[3]*v[1])%p];

\\ Verify: a^3 = [36,0,0]
a_  = [0,1,0]; a2_ = [0,0,1];
print("a^3 = ", fmul(a_, a2_), " (expect [36,0,0])");
\\ beta = 2*alpha: beta^3 = 8*36 = 30 = -13 mod 43
b_ = fscl(2, a_);
print("(2a)^3 = ", fmul(b_,fmul(b_,b_)), " (expect [30,0,0])");

\\ zeta3 = 6 in F_43: 6^2+6+1 = 43 = 0 mod 43
z3 = 6; z3sq = (z3*z3)%p;
print("z3^2 = ", z3sq, "  z3^3 mod p = ", (z3^3)%p);

\\ === Sigma_0 matching: (alpha,beta),(alpha*z3,beta*z3),(alpha*z3^2,beta*z3^2) ===
\\ s0 = alpha+beta = 3*alpha,  q0 = alpha*beta = 2*alpha^2
s0 = fscl(3, a_);
q0 = fscl(2, a2_);

\\ G_i = x^2 - z3^{i-1}*s0*x + z3^{2(i-1)}*q0  (stored as [const, x_coeff, x^2_coeff])
G1c = q0;                         G1x = fneg(s0);                          G1x2 = fcon(1);
G2c = fscl(z3sq,q0);              G2x = fneg(fscl(z3,s0));                 G2x2 = fcon(1);
G3c = fscl(z3,q0);                G3x = fneg(fscl(z3sq,s0));               G3x2 = fcon(1);

\\ Derivatives: G_i' = 2x + G_ix  (stored as [const, x_coeff])
G1px = fcon(2); G1pc = G1x;
G2px = fcon(2); G2pc = G2x;
G3px = fcon(2); G3pc = G3x;

\\ Delta = det of [G1x2,G1x,G1c; G2x2,G2x,G2c; G3x2,G3x,G3c]
\\ = G1x2*(G2x*G3c-G3x*G2c) - G1x*(G2x2*G3c-G3x2*G2c) + G1c*(G2x2*G3x-G3x2*G2x)
D0 = fadd(fadd(fmul(G1x2,fadd(fmul(G2x,G3c),fneg(fmul(G3x,G2c)))),fneg(fmul(G1x,fadd(fmul(G2x2,G3c),fneg(fmul(G3x2,G2c)))))),fmul(G1c,fadd(fmul(G2x2,G3x),fneg(fmul(G3x2,G2x)))));
print("Delta_0 = ", D0, " (expect scalar [39,0,0]: 3*(6-36)*3*2*36 mod 43)");
print("  verify: 3*(-30)*6*36 mod 43 = ", (3*((z3-z3sq)%p)*3*2*36)%p);

\\ Delta_0 should be [39,0,0]; invert mod 43
D0inv = fcon(lift(Mod(39,p)^(-1)));
print("Delta_0_inv = ", D0inv, " (= 32 mod 43 since 39*32=1248=29*43+1)");

\\ H1 = (G2'*G3 - G2*G3') / Delta_0  (quadratic in x)
\\ G2'*G3 coeffs:
\\   x^2: G2px*G3x2 + G2pc*G3x2 ... actually:
\\   G2' = G2px*x + G2pc; G3 = G3x2*x^2 + G3x*x + G3c
\\   x^3: G2px*G3x2  x^2: G2px*G3x+G2pc*G3x2  x^1: G2px*G3c+G2pc*G3x  x^0: G2pc*G3c
\\   G2*G3' = (G3px*x+G3pc): x^3: G2x2*G3px  x^2: G2x2*G3pc+G2x*G3px  x^1: G2x*G3pc+G2c*G3px  x^0: G2c*G3pc
\\ H1 num x^3 = G2px*G3x2 - G2x2*G3px = [2]-[1]=0 ✓ (both = [2,0,0])
H1n2 = fadd(fadd(fmul(G2px,G3x),fmul(G2pc,G3x2)),fneg(fadd(fmul(G2x2,G3pc),fmul(G2x,G3px))));
H1n1 = fadd(fadd(fmul(G2px,G3c),fmul(G2pc,G3x)),fneg(fadd(fmul(G2x,G3pc),fmul(G2c,G3px))));
H1n0 = fadd(fmul(G2pc,G3c),fneg(fmul(G2c,G3pc)));

H1x2 = fmul(H1n2,D0inv); H1x1 = fmul(H1n1,D0inv); H1x0 = fmul(H1n0,D0inv);
print("H1: x^2=",H1x2,"  x^1=",H1x1,"  x^0=",H1x0);

\\ H2 = (G3'*G1 - G3*G1') / Delta_0
H2n2 = fadd(fadd(fmul(G3px,G1x),fmul(G3pc,G1x2)),fneg(fadd(fmul(G3x2,G1pc),fmul(G3x,G1px))));
H2n1 = fadd(fadd(fmul(G3px,G1c),fmul(G3pc,G1x)),fneg(fadd(fmul(G3x,G1pc),fmul(G3c,G1px))));
H2n0 = fadd(fmul(G3pc,G1c),fneg(fmul(G3c,G1pc)));
H2x2 = fmul(H2n2,D0inv); H2x1 = fmul(H2n1,D0inv); H2x0 = fmul(H2n0,D0inv);
print("H2: x^2=",H2x2,"  x^1=",H2x1,"  x^0=",H2x0);

\\ H3 = (G1'*G2 - G1*G2') / Delta_0
H3n2 = fadd(fadd(fmul(G1px,G2x),fmul(G1pc,G2x2)),fneg(fadd(fmul(G1x2,G2pc),fmul(G1x,G2px))));
H3n1 = fadd(fadd(fmul(G1px,G2c),fmul(G1pc,G2x)),fneg(fadd(fmul(G1x,G2pc),fmul(G1c,G2px))));
H3n0 = fadd(fmul(G1pc,G2c),fneg(fmul(G1c,G2pc)));
H3x2 = fmul(H3n2,D0inv); H3x1 = fmul(H3n1,D0inv); H3x0 = fmul(H3n0,D0inv);
print("H3: x^2=",H3x2,"  x^1=",H3x1,"  x^0=",H3x0);
print("");

\\ Multiply H1*H2 (degree 4)
\\ [c0,c1,c2] * [d0,d1,d2] = c0*d0 + (c0*d1+c1*d0)x + (c0*d2+c1*d1+c2*d0)x^2 + ...
P0 = fmul(H1x0,H2x0);
P1 = fadd(fmul(H1x0,H2x1),fmul(H1x1,H2x0));
P2 = fadd(fadd(fmul(H1x0,H2x2),fmul(H1x1,H2x1)),fmul(H1x2,H2x0));
P3 = fadd(fmul(H1x1,H2x2),fmul(H1x2,H2x1));
P4 = fmul(H1x2,H2x2);
print("H1*H2 coefficients: x^0=",P0,"  x^1=",P1,"  x^2=",P2,"  x^3=",P3,"  x^4=",P4);

\\ Multiply (H1*H2) * H3 (degree 6)
Q0 = fmul(P0,H3x0);
Q1 = fadd(fmul(P0,H3x1),fmul(P1,H3x0));
Q2 = fadd(fadd(fmul(P0,H3x2),fmul(P1,H3x1)),fmul(P2,H3x0));
Q3 = fadd(fadd(fmul(P1,H3x2),fmul(P2,H3x1)),fmul(P3,H3x0));
Q4 = fadd(fadd(fmul(P2,H3x2),fmul(P3,H3x1)),fmul(P4,H3x0));
Q5 = fadd(fmul(P3,H3x2),fmul(P4,H3x1));
Q6 = fmul(P4,H3x2);

print("");
print("=== H1*H2*H3 (Richelot image sextic) over F_{43^3}: ===");
print("x^0: ",Q0); print("x^1: ",Q1); print("x^2: ",Q2); print("x^3: ",Q3);
print("x^4: ",Q4); print("x^5: ",Q5); print("x^6: ",Q6);
print("");

\\ Check all coefficients are in F_43 (alpha-components = 0)
ok = 1;
if (Q0[2]!=0||Q0[3]!=0, print("ERROR: Q0 not in F_43"); ok=0);
if (Q1[2]!=0||Q1[3]!=0, print("ERROR: Q1 not in F_43"); ok=0);
if (Q2[2]!=0||Q2[3]!=0, print("ERROR: Q2 not in F_43"); ok=0);
if (Q3[2]!=0||Q3[3]!=0, print("ERROR: Q3 not in F_43"); ok=0);
if (Q4[2]!=0||Q4[3]!=0, print("ERROR: Q4 not in F_43"); ok=0);
if (Q5[2]!=0||Q5[3]!=0, print("ERROR: Q5 not in F_43"); ok=0);
if (Q6[2]!=0||Q6[3]!=0, print("ERROR: Q6 not in F_43"); ok=0);
if (ok, print("All coeffs in F_43 ✓"));

\\ Extract scalars
r0=Q0[1]; r1=Q1[1]; r2=Q2[1]; r3=Q3[1]; r4=Q4[1]; r5=Q5[1]; r6=Q6[1];
print("F_43 sextic: ",r6,"x^6 + ",r5,"x^5 + ",r4,"x^4 + ",r3,"x^3 + ",r2,"x^2 + ",r1,"x + ",r0);
print("Z/3Z check: x^5=",r5,"  x^4=",r4,"  x^2=",r2,"  x^1=",r1,"  (all should be 0)");
print("");

\\ Normalize: monic form y^2 = x^6 + a*x^3 + b
lc = r6; if (lc!=0, lcinv=lift(Mod(lc,p)^(-1)); aa=(r3*lcinv)%p; bb=(r0*lcinv)%p; print("Monic: a=",aa,"  b=",bb), print("lc=0: degenerate Richelot"));
print("");

\\ Verify char poly
print("=== Verify Frobenius char poly ===");
hsex = Mod(r6,p)*x^6+Mod(r5,p)*x^5+Mod(r4,p)*x^4+Mod(r3,p)*x^3+Mod(r2,p)*x^2+Mod(r1,p)*x+Mod(r0,p);
if (r6!=0, cp=hyperellcharpoly(hsex); print("Char poly: ",cp); nj=subst(cp,variable(cp),1); print("#Jac=",nj,"  target=1767  match=",nj==1767), print("SKIP: degenerate (lc=0)"));
print("");

\\ Compare against 7 canonical classes
print("=== Match against 7 Z/3Z canonical classes ===");
ca=[1,2,3,4,6,8,16]; cb=[2,8,42,32,39,42,39];
forstep (i=1,7,1, if(aa==ca[i]&&bb==cb[i], print("MATCH class ",i,": a=",ca[i]," b=",cb[i]," ***"), print("no match class ",i,": a=",ca[i]," b=",cb[i])));
print("");

\\ === Try sigma_1: pair (alpha, z3*beta)=>(alpha, 2*z3*alpha), (z3*alpha,z3^2*beta), (z3^2*alpha,beta) ===
\\ s1 = alpha + z3*beta = alpha + 2*z3*alpha = (1+2*z3)*alpha
\\ q1 = alpha * z3*beta = z3 * 2 * alpha^2
s1c1=(1+2*z3)%p; q1c2=(2*z3)%p;
s1 = fscl(s1c1, a_); q1 = fscl(q1c2, a2_);
print("sigma_1: s1=",s1c1,"*alpha  q1=",q1c2,"*alpha^2");

G1s1c=q1; G1s1x=fneg(s1); G1s1x2=fcon(1);
G2s1c=fscl(z3sq,q1); G2s1x=fneg(fscl(z3,s1)); G2s1x2=fcon(1);
G3s1c=fscl(z3,q1); G3s1x=fneg(fscl(z3sq,s1)); G3s1x2=fcon(1);

G1s1px=fcon(2); G1s1pc=G1s1x;
G2s1px=fcon(2); G2s1pc=G2s1x;
G3s1px=fcon(2); G3s1pc=G3s1x;

D1=fadd(fadd(fmul(G1s1x2,fadd(fmul(G2s1x,G3s1c),fneg(fmul(G3s1x,G2s1c)))),fneg(fmul(G1s1x,fadd(fmul(G2s1x2,G3s1c),fneg(fmul(G3s1x2,G2s1c)))))),fmul(G1s1c,fadd(fmul(G2s1x2,G3s1x),fneg(fmul(G3s1x2,G2s1x)))));
print("Delta_1 = ",D1);

H1s1n2=fadd(fadd(fmul(G2s1px,G3s1x),fmul(G2s1pc,G3s1x2)),fneg(fadd(fmul(G2s1x2,G3s1pc),fmul(G2s1x,G3s1px))));
H1s1n1=fadd(fadd(fmul(G2s1px,G3s1c),fmul(G2s1pc,G3s1x)),fneg(fadd(fmul(G2s1x,G3s1pc),fmul(G2s1c,G3s1px))));
H1s1n0=fadd(fmul(G2s1pc,G3s1c),fneg(fmul(G2s1c,G3s1pc)));

H2s1n2=fadd(fadd(fmul(G3s1px,G1s1x),fmul(G3s1pc,G1s1x2)),fneg(fadd(fmul(G3s1x2,G1s1pc),fmul(G3s1x,G1s1px))));
H2s1n1=fadd(fadd(fmul(G3s1px,G1s1c),fmul(G3s1pc,G1s1x)),fneg(fadd(fmul(G3s1x,G1s1pc),fmul(G3s1c,G1s1px))));
H2s1n0=fadd(fmul(G3s1pc,G1s1c),fneg(fmul(G3s1c,G1s1pc)));

H3s1n2=fadd(fadd(fmul(G1s1px,G2s1x),fmul(G1s1pc,G2s1x2)),fneg(fadd(fmul(G1s1x2,G2s1pc),fmul(G1s1x,G2s1px))));
H3s1n1=fadd(fadd(fmul(G1s1px,G2s1c),fmul(G1s1pc,G2s1x)),fneg(fadd(fmul(G1s1x,G2s1pc),fmul(G1s1c,G2s1px))));
H3s1n0=fadd(fmul(G1s1pc,G2s1c),fneg(fmul(G1s1c,G2s1pc)));

\\ Check if Delta_1 is scalar; invert
if (D1[2]==0&&D1[3]==0, D1inv=fcon(lift(Mod(D1[1],p)^(-1))),
  print("Delta_1 not scalar: ",D1); D1inv=fcon(0));

H1s1x2=fmul(H1s1n2,D1inv); H1s1x1=fmul(H1s1n1,D1inv); H1s1x0=fmul(H1s1n0,D1inv);
H2s1x2=fmul(H2s1n2,D1inv); H2s1x1=fmul(H2s1n1,D1inv); H2s1x0=fmul(H2s1n0,D1inv);
H3s1x2=fmul(H3s1n2,D1inv); H3s1x1=fmul(H3s1n1,D1inv); H3s1x0=fmul(H3s1n0,D1inv);

\\ Product sigma_1
P1_0=fmul(H1s1x0,H2s1x0); P1_1=fadd(fmul(H1s1x0,H2s1x1),fmul(H1s1x1,H2s1x0));
P1_2=fadd(fadd(fmul(H1s1x0,H2s1x2),fmul(H1s1x1,H2s1x1)),fmul(H1s1x2,H2s1x0));
P1_3=fadd(fmul(H1s1x1,H2s1x2),fmul(H1s1x2,H2s1x1)); P1_4=fmul(H1s1x2,H2s1x2);

Q1_0=fmul(P1_0,H3s1x0); Q1_1=fadd(fmul(P1_0,H3s1x1),fmul(P1_1,H3s1x0));
Q1_2=fadd(fadd(fmul(P1_0,H3s1x2),fmul(P1_1,H3s1x1)),fmul(P1_2,H3s1x0));
Q1_3=fadd(fadd(fmul(P1_1,H3s1x2),fmul(P1_2,H3s1x1)),fmul(P1_3,H3s1x0));
Q1_4=fadd(fadd(fmul(P1_2,H3s1x2),fmul(P1_3,H3s1x1)),fmul(P1_4,H3s1x0));
Q1_5=fadd(fmul(P1_3,H3s1x2),fmul(P1_4,H3s1x1)); Q1_6=fmul(P1_4,H3s1x2);

print("sigma_1 result (x^0..x^6):");
print(Q1_0,Q1_1,Q1_2,Q1_3,Q1_4,Q1_5,Q1_6);
ok1=1;
if(Q1_0[2]!=0||Q1_0[3]!=0,ok1=0); if(Q1_3[2]!=0||Q1_3[3]!=0,ok1=0);
if(Q1_6[2]!=0||Q1_6[3]!=0,ok1=0);
if (ok1, lc1=Q1_6[1]; if(lc1!=0, lc1inv=lift(Mod(lc1,p)^(-1)); aa1=(Q1_3[1]*lc1inv)%p; bb1=(Q1_0[1]*lc1inv)%p; print("sigma_1 monic: a=",aa1,"  b=",bb1), print("sigma_1: lc=0")), print("sigma_1: non-F_43 coeffs"));
print("");
print("=== Done ===");

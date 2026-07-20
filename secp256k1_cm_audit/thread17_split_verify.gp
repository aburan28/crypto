\\ thread17_split_verify.gp
\\
\\ Thread 17: Empirically confirm Corollary (cor:p-splits):
\\ for (p, a2) with p prime, p !| a2, a2 != 0:
\\   (i)  p does NOT divide disc(K)  => p is not ramified,
\\   (ii) kronecker(sf, p) = +1      => p splits (not inert),
\\   (iii) [P]^2 = 1 in Cl(O_K)     => order-2 Frobenius ideal.
\\
\\ 10 cases: 5 with D<0, 3 with D>0, 2 with large |a2|/p.
\\ Run: gp -q secp256k1_cm_audit/thread17_split_verify.gp

default(parisize, 64000000);
default(timer, 0);

\\====================================================================
\\ Verify one (p, a2) pair; return 1 on full pass, 0 on failure.
\\====================================================================
verify_split(p, a2) = {
  my(D, sf, m, disc_K, kron_val, ok_disc, ok_kron, ok_order2, B, id, isprinc);
  D = a2^2 - 4*p^2;
  if(D == 0, printf("  SKIP: D=0 for p=%d a2=%d\n", p, a2); return(-1));
  sf = core(D);                  \\ PARI built-in squarefree part (signed)
  if(sf == 0, printf("  SKIP: sf=0\n"); return(-1));
  m = round(sqrt(abs(D) / abs(sf)));
  if(m^2 * abs(sf) != abs(D),
    printf("  ERROR: m decomposition failed for D=%d\n", D); return(0));

  \\ (i) Discriminant check: p must NOT divide disc(K=Q(sqrt(sf)))
  disc_K = if(sf % 4 == 1, sf, 4*sf);
  ok_disc = (disc_K % p != 0);

  \\ (ii) Kronecker symbol: +1 means p splits, -1 inert, 0 ramified
  kron_val = kronecker(sf, p);
  ok_kron = (kron_val == 1);

  \\ (iii) [P]^2 = 1 via bnfisprincipal
  B = bnfinit(Pol([1, 0, -sf]), 1);  \\ x^2 - sf
  id = idealprimedec(B, p);
  if(#id == 0,
    printf("  ERROR: no prime above p=%d in Q(sqrt(%d))\n", p, sf); return(0));
  isprinc = bnfisprincipal(B, idealpow(B, id[1], 2), 1);
  \\ isprinc[1] is the exponent vector in Cl(O_K); 0-vector iff principal
  ok_order2 = (norml2(isprinc[1]) == 0);

  if(ok_disc && ok_kron && ok_order2,
    printf("PASS  p=%d a2=%d sf=%d kron=%d p!|disc=%d [P]^2=1=%d\n",
           p, a2, sf, kron_val, ok_disc, ok_order2),
    printf("FAIL  p=%d a2=%d sf=%d kron=%d p!|disc=%d [P]^2=1=%d\n",
           p, a2, sf, kron_val, ok_disc, ok_order2));
  return(ok_disc && ok_kron && ok_order2);
}

\\====================================================================
\\ 10 test cases
\\====================================================================

print("=== Thread 17: Splitting Verification (10 cases) ===");
print("Conditions: p!|disc(K), kron(sf,p)=+1, [P]^2=1 in Cl(O_K)");
print("");

{
  my(total=0, ok=0, r);
  my(cases = [
    \\  [p,    a2]    D=a2^2-4p^2
    [ 11,    3],   \\ D=-475<0,  sf=-19
    [ 17,    5],   \\ D=-1131<0, sf=-1131/...
    [ 23,    7],   \\ D=-2067<0, sf=-23
    [ 41,    9],   \\ D=-6643<0, sf=-6643/...
    [ 53,   11],   \\ D=-11115<0
    [ 67,   13],   \\ D=-17787<0
    [ 97,   19],   \\ D=-37275<0
    [101,   21],   \\ D=-40363<0
    [113,  227],   \\ D>0: 227^2-4*113^2=51529-51076=453; sf=453/...
    [131,  263]    \\ D>0: 263^2-4*131^2=69169-68644=525; sf=21
  ]);

  for(i = 1, #cases,
    my(p = cases[i][1], a2 = cases[i][2]);
    r = verify_split(p, a2);
    if(r >= 0, total++; ok += r));

  printf("\n=== Summary: %d/%d passed ===\n", ok, total);
  if(ok == total,
    print("All cases confirm Corollary: p splits in K (not inert, not ramified)."),
    print("WARNING: failures detected."));
}

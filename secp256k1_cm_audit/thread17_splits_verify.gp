\\ Thread 17: Verify Proposition prop:order2-frob
\\ For 10 new (p, a2) pairs, confirm:
\\   (i)  p splits in Q(sqrt(sf(a2^2-4p^2))): Kronecker(sf,p)=+1
\\   (ii) [P]^2 = 1 in Cl(O_K)

default(parisize, 64000000);
default(timer, 0);

sf_part(n) = if(n == 0, 0, sign(n) * core(abs(n)));

check_one(p, a2) = {
  my(D, sf, K, Pp, P, P2, v, ok, split_ok);
  D  = a2^2 - 4*p^2;
  sf = sf_part(D);
  K  = bnfinit(x^2 - sf, 1);
  Pp = idealprimedec(K, p);
  split_ok = (#Pp == 2);
  if(!split_ok,
    printf("p=%3d a2=%2d sf=%4d NOT SPLIT (#Pp=%d) -- UNEXPECTED\n", p, a2, sf, #Pp);
    return(0));
  P  = Pp[1];
  P2 = idealpow(K, P, 2);
  v  = bnfisprincipal(K, P2);
  ok = (v[1] == 0*v[1]);
  printf("p=%3d a2=%2d sf=%5d h=%d split:YES [P]^2=1:%s\n",
         p, a2, sf, K.clgp.no, if(ok,"YES","FAIL(!)"));
  ok
};

{
  my(ps, a2s, pass, fail, ok, n);
  \\ 10 new test pairs not in Threads 15-16
  ps  = [101, 103, 107, 109, 113, 127, 131, 137, 139, 149];
  a2s = [  7,  11,  13,  17,  19,  23,  29,  31,  37,  41];
  pass = 0; fail = 0;
  n = length(ps);
  for(i=1, n,
    ok = check_one(ps[i], a2s[i]);
    if(ok, pass++, fail++));
  printf("\n=== %d/%d passed (split + [P]^2=1) ===\n", pass, n);
}

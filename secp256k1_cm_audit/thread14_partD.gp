\\ thread14_partD.gp — Frobenius order in CM fields for norm-form primes
default(parisize, 64000000);
default(timer, 0);

get_order(ce, cyc) = {
  my(n, ord, tmp, all_zero, i);
  n = #cyc;
  tmp = Vec(ce);
  if(#tmp != n, return(-1));
  ord = 1;
  while(ord <= 1000,
    all_zero = 1;
    for(i = 1, n, if(tmp[i] != 0, all_zero = 0; break));
    if(all_zero, return(ord));
    ord++;
    for(i = 1, n, tmp[i] = Mod(tmp[i] + lift(ce[i]), cyc[i]))
  );
  -1
};

print_case(k, p, sf4) = {
  my(Kp, h, cyc, Pp, res, ce, ord);
  Kp = bnfinit(x^2 - sf4, 1);
  h = Kp.clgp.no;
  cyc = Kp.clgp.cyc;
  Pp = idealprimedec(Kp, p);
  res = bnfisprincipal(Kp, Pp[1]);
  ce = res[1];
  ord = get_order(ce, cyc);
  printf("  sf=%-9d k=%-4d p=%-8d h=%-4d cyc=%Ps class_exp=%Ps ord=%d\n", sf4, k, p, h, cyc, ce, ord)
};

print("CM-73 reference:");
print_case(1, 19, -219);
print();
print("Non-CM-73 norm-form primes (small |sf| first):");
print_case(21, 349, -939);
print_case(105, 8287, -1731);
print_case(131, 12889, -3);
print_case(25, 487, -3819);
print_case(35, 937, -5619);
print_case(31, 739, -8643);
print_case(119, 10639, -32187);
print_case(65, 3187, -32619);
print_case(109, 8929, -35859);
print_case(99, 7369, -43059);
print_case(41, 1279, -14619);
print_case(55, 2287, -16419);
print_case(85, 5437, -61995);
print_case(91, 6229, -71499);
print_case(101, 7669, -87267);
print();
print("DONE.");

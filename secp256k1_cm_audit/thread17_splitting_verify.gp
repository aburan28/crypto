/* thread17_splitting_verify.gp
 * Thread 17: verify p ALWAYS SPLITS in Q(sqrt(sf)) for biquadratic Weil polys.
 */

check_one(p, a2) = {
  my(D, sf, m, disc_K, leg, ramified, inert, splits);
  D = a2^2 - 4*p^2;
  if (D == 0, error("D=0"));
  sf = core(D);
  m = sqrtint(abs(D \ sf));
  if (sf * m^2 != D, error("sf*m^2 != D"));
  if (sf % 4 == 1, disc_K = sf, disc_K = 4*sf);
  leg = kronecker(sf, p);
  ramified = ((disc_K % p) == 0);
  inert    = (!ramified && leg == -1);
  splits   = (!ramified && leg == 1);
  printf("  p=%d a2=%d sf=%d disc_K=%d kron=%d ram=%d inert=%d split=%d\n", p, a2, sf, disc_K, leg, ramified, inert, splits);
  if (ramified, error("UNEXPECTED RAMIFICATION"));
  if (inert,    error("UNEXPECTED INERTNESS"));
  if (!splits,  error("NOT SPLIT"));
  1
}

print("=== Thread 17: p always splits when p!|a2 ===");

check_one(101, 14);
check_one(103, -20);
check_one(997, 44);
check_one(1009, -100);
check_one(9973, 200);
check_one(99991, 500);
check_one(999983, 1000);
check_one(1299709, -888);
check_one(15485863, 2000);
check_one(32452843, -3600);

print("10/10 splits confirmed");

print("=== Ramification impossibility sweep (p<=200, |a2|<=50) ===");
rcount = 0;
forprime(pp = 3, 200,
  for(aa2 = -50, 50,
    if (aa2 == 0 || aa2 % pp == 0, next());
    my(D2 = aa2^2 - 4*pp^2);
    if (D2 == 0, next());
    my(sf2 = core(D2));
    my(dk2 = if(sf2 % 4 == 1, sf2, 4*sf2));
    if (dk2 % pp == 0, printf("  RAMIFIED: p=%d a2=%d\n", pp, aa2); rcount++)));
if (rcount == 0, print("  None found (ramification impossible, as proved)"), printf("  WARNING: %d found\n", rcount));

print("=== Thread 17 PASSED ===");

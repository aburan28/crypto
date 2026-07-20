/* Ramification impossibility sweep for Thread 17 */
ram_sweep() = {
  my(rcount = 0, D2, sf2, dk2);
  forprime(pp = 3, 200,
    for(aa2 = -50, 50,
      if (aa2 == 0 || aa2 % pp == 0, next());
      D2 = aa2^2 - 4*pp^2;
      if (D2 == 0, next());
      sf2 = core(D2);
      if (sf2 % 4 == 1, dk2 = sf2, dk2 = 4*sf2);
      if (dk2 % pp == 0, printf("RAMIFIED: p=%d a2=%d sf=%d disc=%d\n",pp,aa2,sf2,dk2); rcount++)));
  if (rcount == 0, print("None found (ramification impossible, as proved)"), printf("WARNING: %d found\n",rcount));
  rcount
}

print("=== Ramification impossibility sweep (p<=200, |a2|<=50) ===");
ram_sweep();
print("=== Sweep PASSED ===");

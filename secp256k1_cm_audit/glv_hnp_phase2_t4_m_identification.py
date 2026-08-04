"""
GLV-HNP Phase 2, Thread 23 (part 3): which m produced the logged Exp T4 table?

The 2026-07-29 run #1 log entry records an "Exp T4" K1-sweep table but does NOT
record the number of signatures m used for it:

    12-bit/2557  lam*=0.340  K1 = 2 3 4 6 8 12 16 24  ->  5 5 5 5 5 4 1 0
    12-bit/2677  lam*=0.070  K1 = 2 3 4 6 8 12 16 24  ->  5 5 5 2 0 0 0 0

Part 1 of Thread 23 assumed m came from the module's HIST tuples (8 and 10) and
consequently failed to reproduce the two margin cells, which part 2 (Exp F)
then mis-read as "the logged table is not reproducible".  This script settles
it by sweeping m and calling the ORIGINAL module's own success_rate().

It also measures how noisy the m-response is at 12 bits, which is the reason
Exp T4b ("more data does not rescue it") and Thread 23 Exp C both read as
flat: at these sizes the response to m is non-monotone at the 5-seed level.

Run: python3 glv_hnp_phase2_t4_m_identification.py
"""

src = open('glv_hnp_phase2_lambda_threshold.py').read()
orig = {}
exec(compile(src.split('# ' + '=' * 75)[0], 'glv_hnp_phase2_lambda_threshold.py',
             'exec'), orig)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
SEEDS = [42, 1234, 9999, 555, 31337]
LOGGED = {2659: [5, 5, 5, 5, 5, 4, 1, 0],
          2647: [5, 5, 5, 2, 0, 0, 0, 0]}

print("=" * 78)
print("Thread 23 part 3 -- which m reproduces the logged Exp T4 table?")
print("=" * 78)
print("ORIGINAL module (glv_hnp_phase2_lambda_threshold.py), plain LLL, 5 seeds.\n")

matches = {}
for tag, p, b, n, lam in [("12-bit/2557", 2557, 2, 2659, 1755),
                          ("12-bit/2677", 2677, 2, 2647, 185)]:
    G = orig['find_generator'](p, b, n)
    cur = (p, b, n, lam, G)
    print(f"  {tag}  n={n}  lam*={orig['lam_star'](lam,n):.3f}")
    print(f"    {'logged':>7}: {LOGGED[n]}")
    ms = []
    for m in range(7, 17):
        row = [orig['success_rate'](cur, m, k, SEEDS)[0] for k in K1_GRID]
        hit = row == LOGGED[n]
        if hit:
            ms.append(m)
        print(f"    {'m=' + str(m):>7}: {row}   {'<== MATCH' if hit else ''}")
    matches[n] = ms
    print()

common = sorted(set(matches[2659]) & set(matches[2647]))
print(f"  m reproducing BOTH logged rows: {common}")
print("  => the 2026-07-29 Exp T4 table was measured at m = "
      f"{common[0] if common else '???'}, a parameter the log omitted.")
print("     Thread 23 part 2 Exp F is superseded: the table IS reproducible.")

print("\n" + "-" * 78)
print("Non-monotonicity of the m-response at 12 bits (the noise floor)")
print("-" * 78)
print("This is why Exp T4b and Thread 23 Exp C both read as 'more data does")
print("not help': at 5 seeds the 12-bit m-response is not monotone.\n")
for tag, p, b, n, lam, k1 in [("12-bit/2557", 2557, 2, 2659, 1755, 16),
                              ("12-bit/2677", 2677, 2, 2647, 185, 6)]:
    G = orig['find_generator'](p, b, n)
    cur = (p, b, n, lam, G)
    row = [orig['success_rate'](cur, m, k1, SEEDS)[0] for m in range(7, 17)]
    print(f"  {tag} K1={k1:>2}  m=7..16 -> {row}")

print("\nDone.")

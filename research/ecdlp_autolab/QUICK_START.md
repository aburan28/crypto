# Quick start: picking this up cold

If you're starting work on this directory without prior context,
read this first. It's a 5-minute overview before diving in.

## What this is

`research/ecdlp_autolab/` documents 60+ cryptanalytic experiments on
the LMFDB elliptic curves used in the AutoLab `ecdlp_index_calculus`
benchmark. Six research programs (v1-v6) have been completed; four
plans (v7-v10) outline future work; v11 is the publication agenda.

## What's already known (no need to re-verify)

1. **The benchmark curves are structurally hard.** Every classical
   attack family (Pohlig-Hellman, MOV, GLV, Smart anomalous, GHS Weil
   restriction, Petit-Quisquater smooth-p-1, Semaev/Diem, structured
   FB GB, Wu's method, quasi-subfield, LLL, p-adic L, isogeny graph,
   ...) has been ruled out for these specific curves. See
   `structural_resistance_proof.md` and `non_naive_attacks_result.md`.

2. **Naive Semaev IC cannot beat Pollard rho on prime fields.**
   Yokoyama et al. 2020 proved this; we empirically confirmed at
   22-28 bits. See `yokoyama_lower_bound.md`.

3. **C+pthreads gives 23-86× speedup** over Sage Python at 81-bit.
   Best is `phase21_rho_v3.c`. See
   `research_program_v6_result.md` for honest budget analysis.

4. **80-bit AutoLab benchmark is NOT feasible** in the 8 CPU-hour
   budget with current C-extension. Gap is 5.5×. Path to close:
   v7 (engineering: Montgomery, AVX, larger partition).

5. **No algorithmic breakthrough achieved.** All session work is
   either:
   - Negative result (most directions)
   - Engineering (C extension, ~5.5× short)
   - Methodological (high-fidelity protocol, scaling validation)

## What to read first

1. `final_synthesis.md` — master synthesis (most important)
2. `PROGRAMS.md` — index of all research programs
3. `research_program_v6_result.md` — most recent honest measurement

## What to read if planning new work

- `research_program_v7.md` + `research_program_v7_execution.md` —
  next concrete engineering effort (~3 weeks)
- `research_program_v8.md` — algorithmic frontiers (months)
- `research_program_v9.md` — long-horizon roadmap (years)
- `research_program_v10.md` — moonshots and catalog of known-bad ideas
- `research_program_v11.md` — publication and dissemination

## What to do if you have a "new" attack idea

1. Check `research_program_v10.md` first. If your idea is M1-M15,
   the entry explains why it (probably) doesn't work.
2. Check the existing `RESEARCH_*.md` at the crypto repo root. They
   cover many related directions.
3. If it's truly new, propose a new phase in v8 (algorithmic) or
   v9 (long-horizon).

## Most common pitfalls

1. **Measuring raw ops/sec only.** The v6 update corrects this
   mistake. Always validate with end-to-end scaling at multiple bit
   sizes.

2. **Confusing the "tested" Pollard rho cost (`O(√n)`) with theoretical
   `O(√n)`.** The constants matter; full rho has overhead beyond
   point arithmetic. See `phase21_validate_scaling.sage`.

3. **Implementing F4/F5 from scratch.** Sage's implementation is
   already optimized; specialized custom code in v6 was slower.
   For F5 specialization (v8 Phase 23), use existing libraries
   (Singular, OSCAR, msolve) as starting points.

4. **Forgetting the Yokoyama bound.** Naive IC provably can't beat
   rho. Don't waste cycles re-confirming. Focus on non-naive.

5. **Underestimating engineering effort for C extensions.** v6
   showed: hand-rolled 128-bit modular arithmetic is 144× slower
   than GMP. Don't reinvent; use established libraries.

## Best reproducible benchmark

```bash
# Build
clang -O3 -I/opt/homebrew/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib \
  -o phase21_rho_v3 phase21_rho_v3.c -lgmp -lpthread

# Run baseline
./phase21_rho_v3 1208925819614629469615699 12 5 \
  942402533238732514353214 807227645196494903606326 \
  500000 4    # 4 threads, 500K steps/thread

# Expected: ~7.4M ops/sec on 4 threads
```

## Most surprising finding

The user's open sub-problem "can the structured-factor-base constraint
be made sparser?" has a *clean answer*: yes structurally (X^m - 1 is
sparser than ∏(X - x_i)), but **no algorithmically** (Sage F4/F5 doesn't
exploit the sparsity, and Yokoyama's lower bound shows even custom GB
can't help). This is the most important methodological finding: clean
mathematical structure ≠ algorithmic advantage.

## Where to look for "what's the next big idea"

- **Most likely incremental win**: v7 Phase 22 (Montgomery + SIMD)
- **Most likely algorithmic win**: v8 Phase 23 (F5 specialized for Semaev)
- **Most likely long-shot win**: v9 Phase 38 (cross-disciplinary,
  random matrix / spectral / information geometry)

## What stays the same forever

The Pollard rho `O(√n)` baseline. Every classical attack must beat
this to be a "breakthrough". As long as that baseline holds (which
Shoup 1997's generic group model says it does), we're working in a
constrained problem space.

## Closing the loop

When you do new work, update:
1. `PROGRAMS.md` — index status
2. The relevant `research_program_v*_result.md`
3. `final_synthesis.md` if findings are major

When in doubt: write a `phase_NN.md` documenting the experiment and
its result, even if negative. Negative results compound; the
resistance map's value is in its completeness.

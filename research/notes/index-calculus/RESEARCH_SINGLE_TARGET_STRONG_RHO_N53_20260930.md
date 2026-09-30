# Is the selected-panel rho at n = 53 a strawman? A strong single-target rho against the direct arm

Follow-on to `RESEARCH_STRONG_RHO_LADDER_20260929.md` (PR #966) and
`RESEARCH_CACHEGRIND_N41_20260929.md` (PR #1004). Written and committed **before any
timing or instruction count of any arm of this note on this host**.

## The claim under test, and why it is suspect

`docs/ic/BOUNDARY_TARGETS.md` (`vs_rho`, `end_to_end_dlp`) carries one surviving online
win for index calculus: the **selected direct-routed four-shard n = 53 panel**, median
direct/rho wall ratio **0.8392** (direct 3.630 s, 9.740 core-s, 1.06 GB RSS; rho 4.308 s,
4.307 core-s), five of five pairs on an AMD EPYC 9V74
(`stage-108-routing-selection-archive-20260913`). Priority #3 of the ledger is to *preserve*
that wall median while cutting the 2.26× core ratio.

That rho is the command built by `scripts/run_koblitz_stage106_shard_routing.py`:
`koblitz_rho_fixture 53 0 signed_frobenius 1 packed <RHO_BATCH_SEED> 476811900269`.
Reading `solve_fixture_packed` in `examples/koblitz_rho_fixture.rs` (not running it):

- it is **not a distinguished-point rho**: `distinguished_bits` is 0, every step is inserted
  into and looked up in a `HashMap`;
- its orbit canonicalization is `raw_canonicalize` on the packed polynomial-basis
  representation, the same family as the PR #955 reference that the ladder showed costs
  about 8× the instructions of the strong rung at the batch cell (371.1 B → 46.4 B, from a
  library multiplier, a normal-coordinate canonicalization and batched inversion together;
  I have not measured how much of that applies to this file);
- it inverts once per step.

The ladder's finding was that the verdict against the compact-orbit IC flipped at the
canonicalization rung. The same suspicion applies to this row: a rho whose single-target
cost is dominated by an O(n) basis scan per step may have been the reference for the last
surviving win. The ledger calls it "matched automorphism-optimized"; "matched" is exactly the
word the ladder retracted for PR #955.

## What I had seen before writing this

The ledger row and the archive's `verification.json` (the numbers above), and the source of
`solve_fixture_packed`. A back-of-envelope: ideal steps
√(π r / 2A) ≈ 6×10⁵ for n = 53 (r ≈ 2^44.26, A = 106), so 4.3 s is about 7 µs per step,
which is large next to the ladder's strong rho (2,370 instructions per step, about 0.5 µs).
I had **not** seen a timing or instruction count of the strong rho at L = 1, of the direct
arm, or of the panel's rho on this host. No n = 53 arm has been run before this note is
committed. **One plumbing smoke test was run at n = 41**, not n = 53, to check the runner's
completion criteria: with the same explicit scalar and `batch_seed` 7, the panel's rho P
took 139,685 steps and 1.39 s of CPU, and the strong rho S 83,392 steps and 0.025 s (one
seed, one run, not part of any registered quantity). It says P is much slower per step than
S at n = 41; I expect the same direction at n = 53 and have not measured it there.

## Arms (all whole-process, this VM, n = 53, a = 0; compile time charged to nobody)

- **D, direct** (the ledger's selected arm): `koblitz_rank_fixture 53 0 1 128 <DIRECT_SEED>
  signed_expanded independent pair_pair_parallel_4096 1 476811900269` with the runner's
  `common_environment()` (`KIC_*` switches, `RAYON_NUM_THREADS=4`), cores 0-3. Same
  completion criteria as the panel's verifier: `relation_rank_summary` status `FULL_RANK`,
  `linear_solution_verified`, `all_relations_group_verified`, factor base with
  `selection_uses_scalar_labels == false`.
- **P, the panel's rho**, unchanged: `koblitz_rho_fixture 53 0 signed_frobenius 1 packed
  <RHO_BATCH_SEED> 476811900269`, one thread on core 0.
- **S, strong single-target rho**: `koblitz_rho_batch_ks_strong 53 0 signed_frobenius 1
  <batch_seed>` with `KIC_RHO_RUNG=3 KIC_RHO_LANES=32 KIC_RHO_DP_BITS=4
  KIC_RHO_EXPLICIT_SCALAR=476811900269`, one thread on core 0. This is the ladder's rung R3
  at one fixture (32 independent walks for the one target sharing a distinguished-point
  table, batched inversion across lanes); its only change for this note is an explicit-scalar
  option (`KIC_RHO_EXPLICIT_SCALAR`), added with a unit test, so it attacks the panel's exact
  target. S gets **one core** while D gets four, which is the hardest setting for S.

Binaries: `koblitz_rank_fixture` and `koblitz_rho_fixture` from the merged tree at the
panel's build command (`cargo build --release --locked --example …`); S from the same tree.
Hashes in `SHA256SUMS_inputs`.

## Measurements and rules (fixed now)

- **Timing blocks.** Nine blocks; each runs D, P and S once, in an order drawn per block from
  `random.Random(20260930)`. S uses `batch_seed` 531310 (frozen here), P its panel seed.
  Recorded: wall, user and system CPU of the child (threads included), `walk_steps` for
  P and S, and the completion criteria. `W_a`, `C_a` are medians over blocks.
- **Walk-seed distribution.** A single-target rho's cost is random over its walk. For the
  same target, P is run under 16 `batch_seed`s (the panel's first) and S under 16
  `batch_seed`s (531310 first), each once, single thread, core 0; the seeds are
  the first 8 bytes, little-endian, of `sha256("SINGLE-TARGET-RHO-N53|{arm}|{i}")`, `i` = 1…15
  after the frozen first. This answers whether the panel's 4.308 s is a typical draw.
- **Instruction counts.** One `valgrind --tool=callgrind` run each of D (four threads), P
  (panel seed) and S (seed 531310), whole process, `Ir`. Descriptive: `Ir_D / Ir_S`,
  `Ir_D / Ir_P`.

**Decision (registered, applied to the nine blocks).** `R2 = W_D / W_S` with a 95 %
percentile bootstrap over blocks (10,000 resamples, `random.Random(20260930)`), resampling
block indices and re-taking the three medians.

- **The selected-panel wall win SURVIVES the strong rho** iff the upper bound of R2 is < 1.
- **It DIES** iff the lower bound of R2 is ≥ 1.
- otherwise **UNRESOLVED**.

**Strawman test.** `Q = median wall of P over its 16 seeds / median wall of S over its 16
seeds`, with a 95 % bootstrap over seeds (separate resampling for P and S). **P is a
strawman** iff the lower bound of `Q` is ≥ 1.5 (the ladder's rule of thumb that anything
under 1.5 is inside uncontrolled build-to-build and walk variation).

**Gates (blocking).** G1: every D run meets the panel's completion criteria and every P and
S run recovers the planted scalar 476811900269 (verified `[d]G = Q`). G2: P and S each have
median `walk_steps` within [0.7, 3] × the ideal √(π r / 2A); a rho outside that band is not
behaving as a rho and the strawman test is withheld. G3 (non-blocking): report `W_D / W_P` on
this host against the panel's 0.839; if it is ≥ 1.0 the panel's headline does not reproduce
on this VM and only the D-versus-S comparison is informative.

**What each outcome licenses.**
- SURVIVES: the online wall crossover of the selected panel stands against a strong rho on
  this host when S is limited to one core; it says nothing about equal cores (the ledger's
  core ratio above 2 already stands) or about the full-cost gate.
- DIES: the selected-panel wall crossover is withdrawn against a strong single-target rho on
  this host; the ledger row and priority #3 are corrected to name S, in the ladder's manner
  (an erratum, measurements retained).
- UNRESOLVED: reported with the interval; no ledger change beyond recording it.

**Inadmissible, stated in advance:** changing S's rung, lanes or DP bits, or D's switches,
after seeing any timing; dropping blocks or seeds; scoring on one block; giving S more cores
after seeing results (equal-core comparisons are a separate, later registration); reading the
panel's EPYC numbers as measured here.

## Scope and honesty limits

- One host class (4 vCPU Xeon here; the panel's EPYC 9V74 numbers are quoted, not
  reproduced). Absolute times differ; only ratios measured here are claimed.
- One target (scalar 476811900269), one n; the walk-seed distribution varies the rho's walk,
  not the target, which is what a single-target rho's cost depends on.
- D is taken as the ledger selected it; I do not tune it. "Fresh build" (compiling the
  producer) is not charged to any arm, so this note does not rescore the ledger's 17.6×
  full-cost ratio.
- S is the strongest single-target rho built in this repository, not the best possible;
  after the ladder's R3 the canonicalization scan is still about half its instructions, so a
  win for D would be a win against a reference that can still be made cheaper.
- Class (`AGENTS.md` §3): **accounting** (reference strength); it moves no boundary target
  except possibly by withdrawing one.

---

<!-- Results are appended below by the commit that follows the runs. Nothing above this line is edited after seeing them. -->

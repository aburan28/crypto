# Sort packed matrix-F5 input rows once by initial density

## Frozen hypothesis and reference

The accepted selective-echelon matrix-F5 fast path returns about
13.7 million n24 degree-4 monomials and spends roughly half its
complete call unpacking them. Choosing the lightest *current* row for
every pivot halved the output term count, but rescanning and
recounting rows made reduction slower
(`research/gf2_minweight_pivots_20260929/RESULT.md`).
Sorting packed rows once by their initial popcount before elimination
may produce lighter pivot rows with only one linear count and one
stable sort. Charge the sort to the complete call and reduction phase.
The candidate may change echelon row order, raw fingerprints, output
term count, and reduction word operations. It must preserve rank,
canonical row space, F5 criterion and row-build counts, and correct
polynomials. The hypothesis is a complete-call improvement large
enough to justify an eligible x86 comparison.

The frozen source reference is main commit
`d38a6fa738d93afb906df70477056498c7299dd6`, with SHA-256
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`,
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

## Frozen work and gates

The new arm sets `KIC_F5_INITIAL_DENSITY_SORT=1`; the prior arm sets
it to `0`. Both use selective echelon output, fused row counting,
direct packed rows, direct scalar unpack, forced AVX2 row XOR where
available, and Gray-code table reuse. Disable deferred-above,
word-batch, AVX-512 unpack, and experimental table builders. Use one
Rayon thread, the seven F5 cases from `examples/f4_f2_bench.rs`,
and seed XORs `0` and `badc0de1` for a local ARM64 screen. The primary
case is `f5_n24_m24_d4`. One release binary and separate arm processes
prevent environment options from crossing arms. The timer includes
criterion, row build, sort, elimination, and full row unpack;
fingerprints and fixture construction remain outside it.

First test that preordering preserves rank and canonical row space on
small, random, sparse, dense and partial-word matrices and that full
F5 reports retain their criterion and build counts. For the local
screen, run one warmup per arm, five prior/prior A/A pairs, then five
alternating prior/new pairs per seed. Preserve every process output,
failure, timeout, and OOM, source/binary hashes, host, phase costs,
signatures, terms and counted work. A local rejection is sufficient
if either primary seed has a complete-call paired median below 1.05,
any smaller case falls below 0.95, or correctness fails. This local
screen cannot establish an x86 speed gain.

If the local gate passes, repeat on Linux x86-64 with AVX2/BMI2, one
pinned allowed CPU, the same five A/A and five alternating pairs on
seeds `0`, `badc0de1`, `5eed2026`, `f5c02a28`. Report exact five-pair
bootstrap 95% intervals, A/A ranges and every smaller case. The
requested further 2× passes only if the frozen complete-call median
and its lower interval bound both reach 2.00, every holdout beats its
A/A maximum, no smaller case falls below its A/A minimum, and all
correctness checks pass. An incremental opt-in remains only if the
frozen median and lower bound exceed 1.05 with the same guards.
Otherwise archive the tested patch and negative receipt and remove
the runtime option. This is a matrix-F5 solver-stage diagnostic, not
an IC online-time or DLP speedup.

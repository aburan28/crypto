# Degree-bounded sparse storage for Boolean polynomial rows

Standalone algebra component · 2026-09-28

**Outcome.** I implemented an on-demand monomial index for squarefree Boolean monomials, packed sparse GF(2) rows, incremental echelon rank, and explicit fill-in limits. It never allocates an array over the full \(2^n\) monomial universe. The implementation is independent of a cryptographic solver.

The earlier 30-variable failure mode called for three full-universe indexing arrays and dense rows of \(2^{30}\) bits. A degree-four bound makes the *potential column set* \(\sum_{j=0}^{4}\binom{30}{j}=31,931\) columns. The index is calculated on demand; even the 31,931-entry lookup table is unnecessary. The bound is an assumption to check against the real algebra workload, not an assertion that its computation stays at degree four.

## Mechanism

Represent a monomial by an integer mask; bit \(i\) records variable \(x_i\). Its degree is `mask.bit_count()`. For the one bits at positions \(i_1<\cdots<i_k\), use a contiguous graded-colex index

\[
\operatorname{index}(M)=\sum_{j=0}^{k-1}\binom nj+\sum_{j=1}^{k}\binom{i_j}{j}.
\]

The combinatorial number system has a unique inverse, implemented by binary search for each \(i_j\). Neither conversion enumerates monomials or allocates the full index space. The index width is a packed 32-bit unsigned integer when the column count is at most \(2^{32}\), and 64 bits up to \(2^{64}\). Configurations above the latter limit fail explicitly.

Rows are sorted, duplicate-free `array('I')` or `array('Q')` values. Duplicate occurrences cancel mod two. XOR is a two-pointer symmetric difference. Echelon reduction keeps one sparse pivot row per leading column; it charges input size, intermediate fill-in, and total stored terms against configurable budgets. A refused row leaves the pivot set unchanged.

The supplied `multiply_row` uses Boolean squarefree multiplication, so \(x_i^2=x_i\). It returns both a degree-bounded output row and the number of discarded terms. Discarding high-degree terms changes the algebraic system in general. A caller must not treat the truncated output as equivalent to the original full polynomial when the dropped count is nonzero. An application needs a justified degree strategy or an explicit incomplete/bounded-search result.

## Results from the actual run

| Variables | Maximum degree | Potential columns | Full-universe dense row | Degree-bounded dense row |
|---:|---:|---:|---:|---:|
| 30 | 4 | 31,931 | 128 MiB | 3,992 bytes |
| 31 | 4 | 36,457 | 256 MiB | 4,558 bytes |

For two **synthetic** 256-row matrices, each input row had 64 terms. Sparse elimination reached rank 256 in both cases. At 30 variables it retained 18,544 terms in pivot rows (74,176 packed bytes) and reached a measured Python heap peak of 122,292 bytes. At 31 variables it retained 18,352 terms (73,408 packed bytes) with a peak of 121,460 bytes. Fill-in raised the longest intermediate row to 312 and 248 terms, respectively. These figures exclude the process interpreter and native allocator overhead; they are `tracemalloc` Python-allocation peaks, not total RSS.

In a separate **storage-only** check at 30 variables, 39,001 rows × 64 unique terms required exactly 9,984,256 bytes of packed index payload. Peak measured Python allocation was 13,621,168 bytes. This run held all rows alive and **did not run elimination**. Its size matches the prior row count, but its rows are synthetic and do not model the sparsity, degree distribution, or fill-in of the blocked solver.

The storage check took 11.45 seconds and the 256-row elimination checks took about 0.15 and 0.14 seconds in one Python 3.12.14 run. These are reproducibility context, not a comparison against a previous solver: no corresponding old implementation was timed, and single-run timing is noisy.

## Verification

- Exhaustively enumerated all degree-bounded monomials at \((n,d)=(6,3)\) and \((8,4)\). The rank map has no gaps or collisions, and every monomial round-trips.
- Checked 10,000 sampled indices for each of \((30,4)\) and \((31,4)\), plus endpoints; all round-tripped. A separate \((90,7)\) case verified the 64-bit index transition without allocating its 8,140,616,224-column space.
- Compared the independent/dependent outcome of each of 180 sparse row insertions against a separate dense-integer GF(2) elimination on a small instance. Both found rank 93 and 87 dependent rows.
- Checked duplicate cancellation, symmetric-difference XOR, Boolean multiplication with a dropped-term count, and three distinct budget failures: oversized input, total stored terms, and elimination fill-in. All passed.

The JSON receipt includes source SHA-256, exact counts, payload sizes and measured memory. Run `python3 test_bounded_sparse_gf2.py` in the extracted archive to repeat the checks.

## Architectural boundary

This fixes **index representation and row storage** for genuinely low-degree, sparse input. It does not yet fix the blocked computation: its actual row polynomials might exceed the chosen degree; elimination might create dense rows, overflow the budgets, or require an algorithm beyond simple GF(2) echelon reduction. In particular, a degree-four cap has not been validated for that workload, and no claim of improved solver throughput or correctness is made.

A safe integration contract for any algebra application would supply its actual monomial masks and rows, record original and discarded degree distributions, expose the row and total-term budgets, and refuse to call a truncated or budget-exceeded run complete. Next, run a **read-only trace** of the blocked workload's monomial degrees and per-row term counts. Only if those fit should this component be tried behind a replaceable storage interface; measure actual fill-in and rank against a small dense reference before relying on larger runs.

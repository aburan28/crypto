# Exact graded reuse inside complete Boolean F5 batches

The native producer now implements the [frozen protocol](PROTOCOL.md). For a
fixed quadratic core, the degree-four projection of each selected Macaulay row
is fixed across affine assignments. `HighCache` compiles the exact word-aligned
high-prefix echelon transform with an identity sidecar. On a later assignment
it recomputes the Boolean F5 criterion, checks the complete selected row-label
sequence and full-column support, and verifies that every prefix word is
unchanged. Those checks make the prefix pivot trace identical. The changing
suffix is multiplied by the retained transform using sparse input terms and
packed transform columns; the inherited M4RI implementation then resumes at
the next word. Any failed guard uses a charged fresh-F5 fallback. The first
all-quadratic assignment often lacks complete column support and uses that
fallback itself, while still supplying a valid fixed high-prefix transform.

Development tests compare the **ordered returned polynomials byte for byte**
against the inherited F5 Echelon path at n=12/16/20/24 and independently
recompute the cached prefix. An ignored n=24 batch-32 test provides local
profiling on public development seeds only. Its timings are unisolated and do
not qualify the registered performance gate. The discovery and holdout seeds
in `protocol.json` remain for the frozen Linux x86-64 campaign.

On a clean committed Linux x86-64 AVX2 checkout, run
`bash research/boolean_f5_graded_batch_20261002/run.sh discovery NEW_DIRECTORY`.
The script uses the repository's existing CPU-isolation controller, retains
source, binary, protocol, resource and output receipts, runs native replay,
and seals successful or failed attempts without overwriting them. A holdout
requires a sealed discovery bundle whose complete-F5-batch gate passed. The
GitHub Actions workflow can dispatch either phase and preserve its bundle.

This is an exact Boolean matrix-F5 batch engineering experiment on generated
public systems. It contains no curve or key interface. A positive stage gate
would not establish natural relation yield, independent rank, full index
calculus cost, or a Pollard-rho crossover; those fields remain null here.

# Fresh photographed-curve volcano campaign

**New run authorized after the historical artifacts became unavailable.** The
old experiment is not recovered or pooled. This replacement uses the same exact
photographed source curve with new securely sampled fixtures and a new frozen
identity. The original complete workload remains required.

## Reconstructed and validated

- Complete selected rational 2-isogeny component: **46 curves, 90 directed maps,
  five levels with populations 1/3/6/12/24**. This is not the whole isogeny class
  over every prime degree. New C labels differ from historical labels.
- Exact source: p=72057980207847001, a=18456521102804593,
  b=19958176276690439; #E=72057979853368804, prime subgroup r=18014494963342201.
- Ten distinct generator points modulo sign and ten distinct nonzero uniform k
  values, shared across generators: 100 problems/curve and **4600 transported
  problem/curve cases**. Linux getrandom and unbiased rejection sampling;
  every raw word and accepted/rejected result is retained.
- PARI/GP independently counted points and group invariant factors on all 46
  curves and verified all 90 rational curve/dual identities symbolically.
  Native C++ independently checked primality, all 4600 scalar/order fixtures,
  9000 directed-map fixture pairs and 900 homomorphism/dual controls.
- Unmodified pinned native C backend tests passed: group 376415 checks, BSGS
  2063 checks and rho 134 checks. CMake and Cargo are unavailable; the library
  was built directly with GCC, and no whole-Rust-library test pass is claimed.

The first full-subgroup rho solve and first ten-target BSGS batch passed
independent scalar checks. Those are initial cases, not the full campaign.
Low-level field instrumentation, full-case solver completion and actual full
IC recovery remain open. No fastest-tool, speedup, hardness or complete-cost
claim is made. Local wall and process-CPU times are descriptive L0 records.

## Recoverable inputs and executable source

[artifacts/PUBLIC_MANIFEST.json](artifacts/PUBLIC_MANIFEST.json) identifies the
[public input/source archive](artifacts/public_inputs.tar.xz). It contains all
4600 answer-free P/Q fixtures, exact curve/map records, public rho walk seeds,
native source snapshots and protocol. Expected/solved scalars and raw getrandom
words are excluded. The complete local archive retains them; publication of
that complete archive requires explicit approval and has not occurred.

The fresh solver-fixture SHA-256 is:

`980d931ce1a804684b160bdab9b213b87688ab9100f0cdf1cd14c390e809acb8`

`sha256sum -c artifacts/public_inputs.tar.xz.sha256` checks the public archive.
Inspect its entries with `tar -tJf`, then extract into a fresh directory.
PUBLIC_MANIFEST.json inside lists exact file hashes and excluded private files.
The public solver input contains no expected scalar. Source is available for
review and reconstruction; the original worker's optional checker-data
preflight and raw-RNG audit cannot be restored from this public bundle alone.
No full-private-artifact restoration is claimed.

PARI runtime Debian URLs and hashes are retained locally for certificate replay.
For rebuilding the native library, use the preserved native snapshot as
upstream: `ln -s source/native_snapshot upstream`, then
`bash source/build_native.sh`. Preserve old source and binary bindings when
rebuilding. No timing or replay record is regenerated merely by extraction.

## Checkpointed execution

The [frozen protocol](PROTOCOL.md) preserves V1 and a pre-execution V2 memory
sizing amendment. BSGS uses m=2^27 and a measured 4-GiB table per generator,
reused across ten targets. Setup/query/verification/destruction costs are kept
separate; any derived or amortized sums retain their own labels. Rho uses the
pinned native default negation-map walk and a secure per-instance seed shared
across curves. Neither solver receives the checker scalar.

The native worker saves every target and all returned statistics/resources,
including failures, as atomic records. Complete cases are resumed only with
the same frozen identity. An interrupted BSGS table is rebuilt and that new
setup is charged as a retained new attempt. Field counters and per-step event
streams are not currently collected and are not imputed as zeros.

`source/launch.sh` performs the frozen checksum preflight and runs two BSGS
workers plus two rho workers. It first spreads coverage across curves, then
continues the full 100-problem-per-curve workload. One invocation per campaign
is required: do not start competing launchers against the same job paths.
Creation of STOP requests a checkpointed stop after the current job unit.
Check actual launcher receipts and process state; this file alone does not
claim the full campaign has been launched or completed.

Reports, per-curve and per-problem summaries, heatmaps and compact raw checkpoint
archives are appended as separate snapshots. No missing curve measurement is
invented. Full IC remains a separate implementation/recovery obligation;
relation-yield/regularity diagnostics cannot discharge it.

Historical recovery metadata remains in
[the earlier handoff](../photo_volcano_recovery_20261017/README.md). Its missing
archive problem remains open for that old evidence, but no longer prevents
this newly authorized replacement run.

# Discovery: the selected F5 row skeleton is stable within each quadratic core

The native source-frozen discovery completed all **72 cells and 1,008 public
generated systems**. Every selected-row bitstream and criterion counter passed
native replay. In all **18 preregistered primary cells** (n=16/20/24, both
affine families, three seeds, batch 32), the largest exact signature class
contained **32/32 assignments**, exceeding the 80% structural screen.
The three seeds at each size produced three distinct full signatures, so
this is reuse *within* one fixed quadratic core, not a universal F5 signature.

| Variables | Systems | Selected labels per system | Pruned labels per system | Koszul / Frobenius prunes | Distinct signatures across the three cores | Primary cells passing |
|---:|---:|---:|---:|---:|---:|---:|
| 12 | 252 | 870 | 78 | 66 / 12 | 3 | bridge only |
| 16 | 252 | 2,056 | 136 | 120 / 16 | 3 | 6/6 |
| 20 | 252 | 4,010 | 210 | 190 / 20 | 3 | 6/6 |
| 24 | 252 | 6,924 | 300 | 276 / 24 | 3 | 6/6 |

These counts are exact row-label and criterion-work diagnostics, not
independent relations, row-space ranks, solve times or complete F5 costs.
The generator quadratic support, generator degrees and multiplier universe
stay fixed. The affine tails vary by independent resampling or one-slot
walks; the criterion's exact prune set is unchanged within each core on this
frozen grid. The aggregate criterion work is 367,416 GF(2) word XORs across
the 1,008 evaluations; individual work can change despite identical
selection bits. A later cache would still pay criterion recomputation unless
it separately proves and prices a sound shortcut.

The unchanged source commit is `7b34fdd5a862a3439ae379db69bdae7877c2fcbe`.
The sealed [discovery bundle](qualified_discovery_01/manifest.json) has
manifest SHA-256 `529255bac0bac72d59d77155fba6f50b476292748431f229eb440f33c0d69564`,
raw SHA-256 `8d72961e39859e8c0d9dd9e132a1df202a83bc464c7457bfcd3b8c1d471e8612`,
and report SHA-256 `27b972f3dba9a0918ec9ba40699abc73526be6fbfc9d5dc8bbb5fda92e93748c`.
All 21 manifest members and a complete native replay passed after copying.
The bundle includes the compiler/host/binary hashes, source and protocol,
full selected-bit hex strings, per-system counters and exact class decisions.

**Decision:** advance to the three untouched holdout seeds under the same
worker and verifier source. This is a structural feasibility result only.
Full F5 timing, output materialization, natural relation yield, IC cost and
rho ratio remain null. No curve or key input is present.

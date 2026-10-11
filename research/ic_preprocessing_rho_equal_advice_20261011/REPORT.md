# Preprocessing rho at equal advice against the IC online class (calculation)

This is lead A of
[`RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md`](../../RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)
§2 and §9 item 1. The ledger lists it as not yet run
(`docs/ic/BOUNDARY_TARGETS.md`, preprocessing-rho control).

**Status: derived, not measured.**
- **What was run.** Nothing new on a curve. Every input is a recorded value, copied with
  its source into [`inputs.json`](inputs.json).
- **How the table was made.**
  [`examples/preprocessing_rho_equal_advice.rs`](../../examples/preprocessing_rho_equal_advice.rs)
  computed [`table.md`](table.md) from `inputs.json` with the note's own model.
- **Why no run.** A real Bernstein–Lange (BL) run would need a curve-based,
  Frobenius-folded implementation first: the repo's `preprocessing_rho.rs` is a `Z_p*` toy.
  It would also need about 10¹¹ precompute steps at n = 73 (below).

## Model (the note's)

- **Folded group size:** `N = r/(2n)`.
- **Advice:** `S` = the IC root-table entries, which the IC builds before any target.
- **BL cost:** online `≈ 2√(N/S)` walk steps; precompute `≈ √(N·S)` steps.
- **The generic preprocessing line:** `S·T² ≳ N` (Corrigan-Gibbs–Kogan).
- **Advice is compared in entries, not bytes.** An IC slot is a 64-bit key plus a 64-bit
  value; a BL entry is an x-coordinate plus a log. The two are within about 2× and were not
  normalised (`docs/ic/PLAN_IC_ACCOUNTING_FIXES_20261007.md` asks for bytes).

## Result

| rung | `log₂ N` | `log₂ S` | BL online steps | IC target probes | IC/BL (ops) | BL online ms (projected at the recorded rho rate, R1–R3) | IC target ms (measured, R1–R3) | `S·T²` above `N` (IC) |
|:--|--:|--:|--:|--:|--:|:--|:--|--:|
| n = 73 | 49.07 | 24.65 | ≈ 9,500 | 457,561 | **≈ 48×** | 41, 147, 147 | 249, 326, 544 | 2^13.2 |
| n = 83 | 45.55 | 24.83 | ≈ 2,600 | 8,845,441 | **≈ 3,400×** | 26, 19, 15 | 7,939, 7,746, 5,623 | 2^25.4 |

1. **The IC online class does not survive the matched control.** At equal advice, BL's online
   stage is about 48× fewer operations at n = 73 and about 3,400× fewer at n = 83.
   - In projected time the margins are about 2–6× at n = 73 and 300–410× at n = 83 (run by run), using
     the paired rho walks' own step rates.
   - IC sits above the generic `S·T² ≈ N` line by 2^13.2 and 2^25.4.
2. **Correction to the note at n = 73.**
   - The note's "IC online probes 5.7e7" and "3,000–6,000×" use the rank stage's mean probes
     per relation (`rank_probes_mean` = 56,896,975 in `ic_summary.jsonl`).
   - The recorded target stage is 457,561 probes in 249 ms. So the margin is about 48×, not
     thousands, and the line is 2^13, not 2^27.
   - At n = 83 the note's figures stand.
3. **Correction to the note's "BL precompute takes seconds to minutes".**
   - BL precompute is `√(N·S)`: about 1.2·10¹¹ steps at n = 73 and 3.9·10¹⁰ at n = 83.
   - At the recorded rho rates of 0.6–2.3·10⁵ steps/s per core, that is days on one core.
   - The IC rank stage is cheaper in operations: 600 × the mean probes, so 3.4·10¹⁰ at
     n = 73 and 1.4·10⁹ at n = 83.
   - For a single target, neither beats plain rho: 3.0·10⁷ ideal steps at n = 73 and
     9.0·10⁶ at n = 83.

**What this changes.**
- The ledger's online-after-precompute rows must be compared with preprocessing rho, not
  plain rho. Against it they lose by one to three orders of magnitude online.
- On total single-target cost both lose to plain rho.
- The ledger is not edited here. Promoting this needs a measured BL run, per the ledger's
  own rules.

## Limits

- The BL online constant (2) and the precompute constant are the note's. Real walks differ by
  small constants.
- A probe and a walk step are different operations: the ops ratio is not a time ratio. The
  projected times use rho rates recorded on a loaded host, and R2/R3 at n = 73 ran about 3.6×
  slower than R1.
- No BL implementation on a curve was run.

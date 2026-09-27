# Narrow arithmetic in the prime-field F4 elimination phases

## Hypothesis

The dense residual matrix and the temporary pair-row accumulator in `f4_fp.rs` hold values below `p < 2^32` in `u64` words. The existing dense row update also branches on each nonzero pivot entry and reduces every product through a 64-bit Barrett quotient. For the benchmark primes `p <= 65521`, a product plus one accumulator value fits in `u32`. An opt-in `u32` accumulator and residual matrix, 32-bit Barrett quotient, and branch-free per-column update may lower both matrix traffic and arithmetic latency. A deferred-reduction variant is a separate follow-on if the narrower branch-free path alone is insufficient; it must be measured as a distinct arm.

## Frozen reference and workloads

- Base revision: `4c3f64b4a3bc0b9b6cf222966c7f4c7426d6dd50` (`origin/main` at branch creation). The original `f4_fp.rs` path remains the reference, selected with `F4_FP_NARROW=0` in the same release binary. `F4_FP_NARROW=1` selects the candidate only for `p <= 65536`; larger primes retain the original path.
- Run `examples/f4_fp_bench.rs` at one repetition with its standard and `large` fixed cases. Primary case: `quad_n8_p65521`, which exercises the residual elimination. Keep every emitted case, including solve cases, for correctness and regression controls.
- Add optional hexadecimal seed XOR to the benchmark's fixed system seed, with zero preserving the original cases. Fresh planted-system holdouts use `00000000badc0de1` and `000000005eed2026`. The exact input generator, case list, seed XORs and source hashes belong in every receipt.
- Use one release binary for both arms on one pinned Linux x86-64 runner, `RAYON_NUM_THREADS=1`, one warmup per mode, five reference/reference A/A pairs, then five alternating-order reference/candidate A/B pairs for each of the three seed workloads. Record all calls, failures and timeouts, process and in-engine F4 wall time, output basis fingerprints, solving fingerprints and report shapes. Match all outputs before calculating ratios. Give each call 180 seconds and stop on any timeout or mismatch.
- Record the CPU model/features, memory, Rust version, source/binary SHA-256 hashes, affinity and load. Treat wall time on a shared virtual runner as noisy. Report A/A range, paired median/minimum and the five-pair exact bootstrap 95% interval. Do not compare absolute milliseconds between different runners.

## Success, stop and accounting

The opt-in path must match every basis and solve fingerprint, number of steps, matrix shape, basis size and verdict. To promote it, the primary `f4_ms` paired median reference/candidate ratio must be at least 1.2 with a bootstrap 95% lower bound above 1; both holdout primary ratios must exceed their A/A noise ranges, and smaller cells must not regress beyond their own A/A noise. If an arm fails or the speed gate is missed, retain the raw result and keep the original path as default.

The measured interval is the complete degree-bounded F4 call, including symbolic preprocessing, both elimination phases and final interreduction. Solve calls are additional diagnostics. This is an internal point-decomposition solver-stage experiment, not an IC candidate or one-target DLP comparison. Online IC time, rho reference, speedup, and scoreboard values remain unknown; no end-to-end gain will be claimed from this stage alone.

## Deferred-reduction follow-on, frozen before its run

Run `36293775047` found exact output and a 1.11× primary gain for mode `1`, below the 1.2× gate. The next arm is `F4_FP_NARROW=2`, compared with the original mode `0` using the same complete standard and large case list, frozen seed, two holdout XORs, one pinned CPU, warmups, five A/A pairs and five A/B pairs. Mode `1` remains in the first receipt as an intermediate control. The follow-on keeps 32-bit stored rows and branch-free inner updates, but accumulates several pivot-row contributions in 64-bit scratch and reduces only when a pivot coefficient is needed or a row is stored. Normalize scratch after at most 256 nonzero pivot updates so the accumulator stays far below `2^64` even at `p=65521`. Apply this to both sparse reducer elimination and residual dense RREF. Row-oriented pivot discovery may change intermediate pivot order, but the unique final RREF, complete F4 basis and all benchmark fingerprints must match mode `0` exactly.

The same 1.2× primary complete-F4 median and 95% lower-bound gate applies to mode `2`; both holdout primary cells must clear their own A/A noise, and smaller cases must not show a material regression. Preserve the first run and all failures. A gain confined to one elimination phase or F4 alone remains a solver-stage diagnostic, not a one-target IC claim.

## Four-thread control frozen after the mode-2 one-thread result

Run `36296331241` cleared the mode-2 one-thread gate with exact outputs. Before default promotion, rerun the same 13 cases and three workloads with `RAYON_NUM_THREADS=4`, pinning four allowed CPUs; keep warmups, five A/A and five alternating A/B pairs. The complete F4 call on the primary case must not regress beyond its own four-thread A/A noise in any workload, smaller cases are regression controls, and every fingerprint and counter must match. A four-thread gain is not required. The runner may differ, so compare paired ratios only within that run. If this control passes, select mode `2` by default for `p <= 65536`, retain `F4_FP_NARROW=0` as reference/rollback and mode `1` as intermediate control, and repeat both thread counts with explicit modes before merging.

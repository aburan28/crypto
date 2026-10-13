# Scope audit of the bounded N83 rho optimization

This is a source and correctness audit, saved before campaign timing. It documents what the new worker can test. It does not report a tuning winner or claim a global optimum. The exact workload is the legacy a=1 curve over GF(2^83) with tail [0,2,4,7], subgroup r=8569786107849059, eigenvalue2998004823950457 and signed-Frobenius action166. It is distinct from the registered tail45 model and degree83/a=0 workload.

## Reference already optimized

The existing strong library already uses distinguished points, the signed-Frobenius quotient in a normal basis, hardware carryless field multiplication when available, and lockstep Montgomery batch inversion. The strong defaults are lanes32/DP4/escape-disabled/cap-factor2000 (src/cryptanalysis/koblitz_strong_rho.rs:440–467). The minimum-reference contract admits this strong control and excludes the older packed control (docs/ic/boundary_targets.json:209–224). The original default-strong study remains unchanged.

Current canonicalization is already more advanced than the older September profile: NormalBasis::least_rotation searches longest cyclic zero runs and compares their candidate rotations (koblitz_strong_rho.rs:639–687); it does not normally scan all83 rotations. The worst case and exact tie handling remain part of its semantics. Linear basis changes use byte tables (:534–559), and canonicalization applies the matching ±lambda^k scalar multiplier (:1056–1090). The September statement that a full linear rotation scan consumed half the instructions is historical evidence and cannot describe today's wide-worker profile without a fresh profile.

All32 jump points are G-only and target independent (:1093–1105); online starts contain Q (:1119–1125). Lanes are independent walks within one operating-system worker, rather than32 CPU threads. Active lanes' denominators are inverted together (:1219–1236). Optional fruitless escape doubles the point and both scalar coefficients while keeping their invariant (:1258–1275); those extra additions/canonicalizations are charged even though walk_steps does not count the extra doubling.

## Sharing IC's field implementation

Both backends execute the same generic StrongRhoG<F,u64> production walk. FoldField only delegates degree, multiplication, squaring and batch inversion to the existing Gf2_128 implementation (examples/n83_rho_optimization_eval.rs:93–114). It changes field engineering without replacing canonicalization, collision recovery, start generation, scalar arithmetic or the point law.

WideGf2 already uses three hardware64-bit carryless products for a128-bit Karatsuba multiplication; reduction folds the sparse tail through XOR/shift iterations and square calls multiply_self (wide_gf2m.rs:138–180). Gf2_128 has runtime hardware dispatch and a tail-fold path. Legacy n83 has tail-degree7, meeting d<=63, 2d<n+2 and n+d<=129; hence its CLMUL fold path applies on the current PCLMULQDQ host (semaev_decomp.rs:845–852,920–930,1051–1061). Its square spreads bits and then reduces (:1063–1078). Giving rho these existing IC primitives removes a preventable implementation asymmetry; it does not guarantee they are faster for every operation mix.

Batch inversion is mathematically faithful: both retain zero entries, multiply all nonzero inputs into one accumulator, invert it, and recover individual inverses with prefix products (WideGf2:214–233; Gf2_128:1115–1133). In the strong loop, batch inversion delegates to the chosen field (:1235–1236). Normal-basis rotation inversion used by standalone starts, exceptional additions and scalar multiplications remains the strong library's inv (:958–977); the batched field inversion itself uses each field's polynomial squaring chain. Faster primitive paths may have different numbers of primitive instructions. Equal StrongRhoCharges means equal reported walk/group bookkeeping, not an equal field-operation or instruction ledger.

Before timing, four new-example tests passed. Exact n83 tests cover102 deterministic field samples including0,1,top-bit and maximal field word; every pairwise product, square and inverse; empty/zero-containing batches with widths0/1/16/32/64/102;32 seeded jump points; selected canonical forms; and one complete separate planted-Q solve. Wide and Fold returned identical recovered scalar, walk_steps, starts, fruitless events, escapes, capped walks, wasted merges, table size, ideal count, action size and StrongRhoCharges. The test Q=12345G and fixed test seeds are correctness fixtures, not training or held-out measurement cases. Campaign same-input kernel cases must independently enforce equivalence too; one unit-test solve is not exhaustive proof of every future trajectory.

## Exact scalar arithmetic and compilation

The subgroup is approximately52.928bits and fits u64. MulMod selects a floating quotient only when m<2^50 (koblitz_strong_rho.rs:699–724). Here it takes exact u128 multiplication/remainder. Enabling that approximate path by increasing the threshold would require a new error bound and correctness implementation; it is excluded. Future exact Barrett/Montgomery scalar reduction or conditional subtraction could be separate implementation candidates, but are not hidden inside this study.

Both backends are monomorphized into the same new release executable from one standard cargo release invocation. The captured commands use --release --locked, default features, one toolchain and one target directory. No target-cpu, LTO, PGO or extra compiler flags were supplied. No Cargo.toml release-profile override or applicable Cargo config was found; relevant flag variables were unset when the evidence was saved. Runtime feature detection selects hardware multiplication for both. The exact captured commands and hashes are in checks/worker_commands.json. Compiler instruction layout/inlining can still affect runtime; kernel comparisons use Wide and Fold in this same executable, not an unexplained ratio of separate old/new builds. The unchanged original worker remains the historical baseline.

## Bounded search omissions

The strict new configuration admits two field backends, lanes16/32/64, DP bits4/6, escape false/true and cap-factor2000. The frozen protocol first compares kernels at the default policy, then searches twelve policies on the chosen backend. This is a staged bounded search.

- Lanes1/8/128 and other widths are omitted. The library permits other positive widths; there is no theorem establishing16/32/64 as the optimum. Width changes inversion amortization, scheduling, active-lane utilization and work abandoned on recovery.
- DP bits0–3,5 and8+ are omitted. Lower DP density reduces table traffic/memory while lengthening walks and collision delay. Abandonment, short cycles and the per-walk length cap8*2^DP change with DP policy. Escape-on may matter more for long walks; a narrow-rung test at DP10 exists (:1645–1683), but is not proof of optimality for this exact legacy curve. Current r53, rather than degree83 alone, sets the scale: ideal work is approximately9,005,148.64 steps, and DP4 ideal-density storage is about562,822 entries. The separately larger a=0 subgroup's memory requirements do not carry over.
- Jump count stays32 (JUMPS at:56–57). Different r-adding maps, jump distributions, hash functions and partition policies are omitted. Selecting a favorable jump or walk seed is prohibited; all declared corpus seeds must contribute.
- Kernel/policy interaction is not exhaustively searched. A backend selected at lanes32/DP4/escape-off might lose at another width/density. The full Cartesian product has24 configurations; the staged search deliberately evaluates18 kernel cases plus108 selected-backend policy cases. Its winner is selected within this procedure, not proven best across every combination.
- The current basis-change tables and canonicalization algorithm are retained. Fully normal-coordinate point storage, other bases, cheaper hash tables, specialized scalar arithmetic, alternative escape logic, equivariant sigma walks, parallel shared-DP implementations, vectorization, GPU implementations and compiler/LTO/PGO exploration are outside this grid. An ECC2K sigma implementation for the a=0 family cannot be silently substituted into the a=1 legacy field/subgroup.
- The fixed cap-factor is not tuned to suppress failures. Library cap returns None without an outcome/charge ledger; missing counters remain unknown. Process timeout, OOM, invalid scalar and mismatched target are retained by the runner and cannot be treated as fast solves.
- Kernel equivalence controls isolate one engineering change. Policy changes can change trajectory, useful collision progress and total step count, so neither fastest seed nor ns/step alone is a valid policy-selection metric.

The protocol selects by mean whole-process cost with a default preference within2%, then freezes one configuration before disjoint held-out targets/seeds. An optimum in that cold-selection metric is not an optimum in prepared-online cost. No benchmark was run for this audit, and no worker/library/frozen-source edit was made during this evidence follow-up.

## Pinned sources and captured validation

Source files below are unchanged from the successful worker validation. SHA-256, rather than a filename or curve degree alone, binds the implemented arithmetic.

| Source | SHA-256 |
| --- | --- |
| examples/n83_rho_optimization_eval.rs | fb14072f4fa5d08e18d8b79130b6d58a220adae0f76369cc9be949691c83eaa0 |
| src/cryptanalysis/koblitz_strong_rho.rs | 3e45be340df4fe5b73998b8f907ad0a1a77d58c06609c73f77f17f2babd90108 |
| src/cryptanalysis/wide_gf2m.rs | 086c786cdb2410c852c4242abc543255776a5ec3fff873ab972d08978907c701 |
| src/cryptanalysis/semaev_decomp.rs | be551d55d574b89bab1e9cd6fe1b1ad298966dd4fcf48a7cd1262d3e5ef4925a |
| Cargo.toml | ae3ac94705e6e3bbda0d67da340800ebad492d2151edb3767e2ab4ae00da31c5 |
| release worker executable | 21827dcc645fd4591b39dd913619f047d267a98d541854bd7c83ee7143d8bbcb |

The worker source remains untracked during preparation; the revision observed when saving this audit was25a56bc8fb94f2ea2ee2d6c790cb0e29b8159119. Source hashes and executable hash are the validation pins. N83_SOURCE_COMMIT was not injected into the original build, so compiled_source_revision is null rather than an invented revision.

Validation output was captured by the execution tool, then transcribed verbatim into checks/worker_tests.log and checks/worker_build.log without rerunning commands. Test result:4 passed,0 failed,9.30 seconds; cached test build1.10 seconds. Release build passed in1.14 seconds. Eleven existing library warnings remain visible. The earlier rustfmt call without RUSTUP_HOME failed because the default home was read-only; checks/worker_format_initial_failure.log preserves that output. With the declared toolchain environment, formatting passed silently. These are correctness/build observations, not benchmark receipts. The root agent owns the separate complete-library-suite record.

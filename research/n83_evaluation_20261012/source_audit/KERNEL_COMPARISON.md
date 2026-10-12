# Stock wide IC and strong rho field kernels

Read-only source audit during the active degree-83 panel. No measured code, frozen files, compiler, tests or benchmarks were changed or executed. Baseline: aburan28/crypto at a8194d616ed8c04489ec46177d369f1b1670360a.

## Finding

**Both stock wide IC and strong rho select hardware PCLMULQDQ on this x86 host.** This is not a native-rho versus portable-IC comparison.

The strong producer reports arithmetic_kernel = pclmulqdq. The stock IC uses the same runtime x86 feature detector and the exact legacy modulus. Therefore its hardware branch is selected under the reported host capabilities. This is a source/dispatch audit, not a new timing measurement or inspection of emitted assembly.

| Operation | Stock wide IC | Strong wide rho |
| --- | --- | --- |
| Field type | semaev_decomp::Gf2_128 | wide_gf2m::WideGf2 |
| Runtime detection | is_x86_feature_detected!("pclmulqdq") | Same detector |
| Unreduced multiplication | Hardware u128 Karatsuba over three 64-bit CLMUL products | Same three-product hardware Karatsuba construction |
| Reduction at legacy n=83 | Short-tail carryless folding | Repeated sparse shifts/XORs until high part vanishes |
| Squaring | Bit spreading followed by carryless folding | mul(a,a), including hardware product and sparse reduction |
| Portable fallback | Software XOR-shift schoolbook, used only without detected hardware; test-only portable constructor | Software 64-bit carryless products when hardware is absent |

Different reduction, squaring, inversion and batching arrangements remain implementation differences. No relative primitive cost has been measured here; identical instruction counts or primitive costs must not be inferred from shared CLMUL use.

## Exact call-path evidence

- examples/koblitz_orbit_dlp_fast.rs:2000–2001 constructs FastBinaryCurve128 and Gf2_128 from the legacy irreducible polynomial.
- src/cryptanalysis/koblitz_fast_arith.rs:441–463 stores Gf2_128 inside FastBinaryCurve128 and calls Gf2_128::new.
- src/cryptanalysis/semaev_decomp.rs:843–852 detects PCLMULQDQ at runtime. The constructor's short-tail conditions at 920–930 select fold_tail when hardware is present.
- Here n=83 and tail degree d=7: d<=63, 2d=14<n+2=85, and n+d=90<=129. The nonzero tail is 1+z²+z⁴+z⁷, represented by 149.
- semaev_decomp.rs:960–990 implements mul_fold/reduce_fold. The public mul at 1051–1062 chooses that path when fold_tail is nonzero. clmul at 1013–1024 independently chooses hardware; its software alternative is not forced by a production option.
- semaev_decomp.rs:1069–1087 squares by bit spreading and folding. The portable constructor at 855–860 is cfg(test).
- src/cryptanalysis/wide_gf2m.rs:73–81 detects PCLMULQDQ; 45–51 invokes _mm_clmulepi64_si128; 127–145 chooses hardware and forms its Karatsuba product. Reduction is 149–169; mul/sqr are 173–181.
- src/cryptanalysis/koblitz_strong_rho.rs:157–173 forwards wide field multiplication, squaring and batched inversion to WideGf2. The exact-point worker constructs WideGf2 through the same factory.

The IC root batching path at example lines 1286–1324 batches denominators and uses gf.mul plus its normal-basis-assisted invert128 routine. It is a different workload from rho's 32-lane group additions and cannot be assigned the same operation cost by counting calls alone.

## Reference gate versus attribution

The current minimum-reference rule in docs/ic/boundary_targets.json, measurement_schema.vs_rho.rho_reference_minimum, requires the strong rho reference or an alternate rho with measured per-step cost no worse on the same target. It explicitly excludes the archived packed/reference rho controls. RESEARCH_FIXED_REFERENCE_RHO_20261001.md:40–42 states that rule.

The paired wall-time contract additionally requires the same public point, matched resource envelope, correct accounting, verified answers and required host/identity/replay receipts. These inspected rules do not require identical field-kernel source between IC and rho. They compare complete implementations of different methods.

Consequently, remaining kernel differences are a limit on attributing a ratio solely to algorithm design. They do not, by themselves, invalidate a correctly qualified whole-implementation comparison. To isolate arithmetic effects, preregister a future matched-backend ablation, retain conversions and all stages, and validate the combined candidate. Such an ablation would be new evidence; do not change the active panel.

Common operation-equivalent S remains unknown without calibrated primitive and phase costs. IC root-index probes and rho walk steps remain distinct native counters.

## Source SHA-256 pins

| Source | SHA-256 |
| --- | --- |
| examples/koblitz_orbit_dlp_fast.rs | 5d6cca46e56a582108fc4bf97e6658a5af4a3a09a11fff3eeec26e4304d58ef1 |
| src/cryptanalysis/semaev_decomp.rs | be551d55d574b89bab1e9cd6fe1b1ad298966dd4fcf48a7cd1262d3e5ef4925a |
| src/cryptanalysis/koblitz_fast_arith.rs | ab156d57e983a8f248bf4fd6c2f6d46d2c38da813a4eadd31e715fc443ad7c1a |
| src/cryptanalysis/wide_gf2m.rs | 086c786cdb2410c852c4242abc543255776a5ec3fff873ab972d08978907c701 |
| src/cryptanalysis/koblitz_strong_rho.rs | 3e45be340df4fe5b73998b8f907ad0a1a77d58c06609c73f77f17f2babd90108 |

The frozen baseline IC executable remains SHA-256 983f28400b7852fc6628ba97d3ea83b2ab9e32d9944658c7eea831ee89a4abf2. Neither it nor the current source kernels were modified by this audit.


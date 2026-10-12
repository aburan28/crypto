# Four-word shared-mask loading: unqualified GPU screen

The single vector4 candidate passed correctness validation but did not meet
the predeclared throughput qualification threshold. The scalar shared-mask
implementation remains the selected RTX PRO 6000 preset. This is a bounded
engineering screen, with no confirmed gain or slowdown.

## Complete-walk measurement

One RTX PRO 6000 Blackwell Server Edition ran both modes at shared1,
compact1, WP2, CLMAD1, batch 16, 256 threads/block, minBlocks 2 and 385,024
workers. Each sample completed **201,863,462,912 scalar updates**, using
1,024 steps and 32 launches. Native code was built with CUDA 13.3.73.

| Screen order | Complete scalar iterations/s | Rate / faster control |
|---|---:|---:|
| Scalar shared control before | 14.640971 B/s | 1.000000 |
| Vector4 | 14.682740 B/s | 1.002853 |
| Scalar shared control after | 14.573526 B/s | 0.995393 |
| Required qualification | 14.714175855 B/s | 1.005000 |

The protocol required the candidate to exceed the faster bracket control
by more than 0.5%. Its observed margin was 0.285288%, and control drift was
−0.460659%. The protocol stopped after five timed rows, including two
excluded warmups. It did not run paired confirmations or timed DP34
collection. The screen does not establish whether a smaller improvement
exists and is not a measurement of a physical throughput ceiling.

Both modes passed all six arithmetic suites, 128 storage cases covering
297,344 records, and the dedicated shared-table probe: 21 scenarios,
21,036 pairs, 114 complete block snapshots and 51,072 words. Full client
replay/restart checks, normalized state comparisons, and 28 checkpoint
children passed before timing. The checkpoint set contained 26 successful
children, two expected worker-geometry rejections and 12 comparisons.
Every timed row completed its exact work budget with zero drops.

Runtime calibration reported 104 registers/thread, zero local bytes,
1,792 function shared bytes and a separate 1,024-byte reservation for both
modes. Source, executable, complete native-code and GPU identity bindings
passed before and after the run. The run made one GPU submission and no
retry: session 15541, app `ap-8wtaL9FddF0MIdQZJffxYc`.

## What was tested

The candidate groups 56 unique mask rows into 14 groups of four, stored as
an aligned `uint4[14][8]` table. It loads a group at first demand and uses
each component at the original 64 mask-consumption positions. Both field
outputs, every arithmetic statement, the field basis and the iteration
function remain unchanged. The table still contains 448 words, and each
helper invocation still reads 224 logical mask bytes.

| Native compiler observation | Scalar shared control | Vector4 |
|---|---:|---:|
| Mask loads | 56 scalar LDS | 14 LDS.128 |
| Complete paired-helper non-NOP instructions | 502 | 460 |
| Selected complete-path visits/scalar | 2,189.75 | 2,147.75 |
| Walk registers | 104 | 104 |
| Stack/local/spill bytes | 0 | 0 |
| Compiled shared extent | 2,816 bytes | 2,816 bytes |

The full compiler review also accounts for the changed cooperative copy,
caller, call sequence and every callee. These observations explain the
candidate's admission to the GPU screen. They do not measure memory
transactions, bank conflicts, execution latency or issue-port utilization;
equal entry register counts do not establish identical live ranges.

## Retained code and evidence

The experimental [candidate patch](candidate.patch) is archived rather than
applied to the production tree. It includes the default-off build option,
generator/header changes, diagnostic markers and the vector-aware device
probe. Applying it to a temporary copy of the seven affected base files
reconstructed every measured candidate file byte for byte. The
[source manifest](source-patch-manifest.json) binds all 160 base and candidate
files.

The [complete GPU result](comparison.json.gz) and
[complete compiler result](compiler.json.gz) use deterministic gzip to
retain their exact original JSON bytes. Their uncompressed SHA256 values are:

- GPU: `ab570300f07c31577132c545ecd47b1181b1254122526e12c14947f3dc5a62df`
- Compiler: `7a8b014039f649ccb16639be7c02fb4a537570f346f4b30449462ec38c60a71a`

The [archive index](index.json) gives compressed and uncompressed hashes,
file sizes and the final disposition. The
[independent GPU review](comparison-review.json),
[compiler review](compiler-independent-review.json),
[call-path supplement](compiler-supplement.json), and
[terminal handoff](terminal-handoff.json) retain the validation details.
The driver, fixed plan, layout gate and final auditor are archived with
their original path assumptions for provenance. The initial reviewer
hex-format error and its narrow correction are retained separately;
the frozen experiment bytes and expected digests were unchanged.

The current 26 B complete-scalar-iterations/s objective remains unachieved.

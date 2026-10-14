# Independent replay and file binding: all 139 degree-83 runs

Observed native validation, 2026-10-12 UTC. See [final_verification.svg](final_verification.svg)
for the retained-stage and verification flow. All 139 returned scalars independently
replayed to the frozen public point. Every planned process and answer was retained;
there were zero invalid answers, zero no-answers and zero failed/capped processes.
The nine reference/Barrett arithmetic pairs matched their exported answer,
jump seed, outcome counters, charges and rho policy. These are recorded-semantic
checks, not a complete intermediate-state proof or independent timing replay.

| Gate | Observed native result |
| --- | --- |
| Frozen process matrix | 10 A/A + 18 arithmetic + 81 policy + 30 holdout = 139 |
| Independent scalar/point replay | 139 of 139 pass |
| Exact eight fresh fixtures and prior exclusion stream | Pass; three training and five holdout, 609 excluded prior signed orbits |
| Reference/Barrett exported semantic pairs | Nine of nine match |
| Invalid answers, missing answers, failures/caps | Zero |
| Worker executable and source identities | One frozen solver executable and one worker source SHA-256 |
| Recorded clean environment and resource envelope | All 139 match CPU3, one thread, 12 GiB address space and 60 s cap |
| Captured stdout versus adjacent producer bytes | All 139 byte-identical |
| File binding checks | 2,439 recorded SHA-256 comparisons; 455 distinct bound files |
| Final retained file snapshot | 1,036 files, including heartbeat and isolation logs |

## Exact identity and independent arithmetic

The model is `icv1-f2m83-t6151469093347-9eba6cd5`: a=b=1,
`z⁸³+z⁷+z⁴+z²+1`, prime subgroup order `8569786107849059`, and generator
`(6614244651398184860292015,8868260410029249391539529)`.
The fresh points and prior signed-orbit exclusion invariant are documented in
[GENERATION_REPORT.md](GENERATION_REPORT.md), with standalone C++ sources and
the frozen generation and fixture-replay receipts.

[check_speed.cpp](check_speed.cpp) uses polynomial shift/XOR multiplication,
Euclidean field inversion and binary-curve affine group arithmetic independently
of the optimized Rust solver. It checks each `dG=Q`, then compares d with the
checker-only sidecar. It binds each point input, config and seed and requires
successful process and producer group-replay flags. All stage and per-fixture
count gates passed. No validation scalar was supplied to a measured child.

The independent replay's own native isolated-run receipt records whole-process
wall `0.886574336 s`, user CPU `0.833018 s`, system CPU `0.051324 s`, and peak
RSS `109188 KiB`, with exit 0 and zero observed contention samples. This work is
outside both timed solver arms and is reported separately; it does not retime
the solvers. The first attempt to use unavailable `/usr/bin/time` was retained
in the study's checks; the successful direct native isolated invocation is in
`records/isolation/verification_attempt2.jsonl`.

## File and metadata provenance

[verify_file_provenance.cpp](verify_file_provenance.cpp) ran after replay. It
compared the actual solver executable, sources named in compiled-source receipts,
original frozen protocol, point files, config files, producer output, stdout and
stderr with retained SHA-256 values. It also checked all 139 observation/execution
bindings, all four pre-launch ordered matrices, actual argv/seed/config against
those matrices, frozen training-selection inputs, CPU-component arithmetic and
RSS-unit conversion. JSON copies of the protocol match structurally; the original
protocol byte hash is bound separately because serializer key ordering can differ.

It separately snapshots every panel, heartbeat and isolation file. The audit
marks which hashes matched a pre-existing recorded binding and which are final
snapshot hashes. A new snapshot hash does not imply that file was sealed before
launch. The successful final build had zero warnings; initial indentation
warnings and their corrected build logs are retained.

The complete answer receipt SHA-256 is
`539751e5dc58c0390b2d36df52e016d674eac51c61d770ffb5ffb6e2597e9988`.
The file-provenance receipt SHA-256 is
`8511411695b7de2a2fd693b64122c42a75dc89bf5b9b694ea8c60af710eb5bd1`.
The frozen protocol SHA-256 is
`f5aa7d7fe2e58ddcfb8177dcc9d9907f4c5aeb859cadea983448b872b1f05d7c`.
Commands and source references are in [README.md](README.md).

## Rejection checks and limits

Controlled negative tests, outside the campaign, rejected a mismatched point-file
hash, a wrong frozen-protocol hash, an existing generation output directory and
a deliberately wrong scalar. An initial wrong-scalar mutation first reached the
producer hash gate; after updating the controlled test clone's hash, the native
group replay rejected d=1 for the frozen target. Test clones alone were modified;
all genuine campaign files and original fixtures remain unchanged. Inputs,
failure receipts and raw logs are in [negative_checks/](negative_checks/).

Verified here: final answer correctness, exact fresh fixture identity and
generation, recorded scalar-kernel equality, complete-count gates and file/metadata
bindings. Not established here: independent reproduction of solver timing,
every intermediate state, a proof-assistant formalization of the signed-orbit
argument, a common calibrated operation unit, independent physical-host conditions
or L2/admission gates. The main report supplies the runtime comparison and its
limited held-out statistical scope. This verification report makes no new speed
or cross-method claim.

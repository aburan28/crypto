# Hosted N83 K600 Gray/hash first paired trial

[GitHub Actions run 38983406418](https://github.com/aburan28/crypto/actions/runs/38983406418)
completed on Linux ARM. The pinned source checkout was
`cb98e85d1c1c1afdd5000bc70eea00969ee86c6c`; the dispatched workflow
head was `8dd40601e80ae8e632229c6b4b633d4f21925513`. The image, cgroup,
source, public corpus, inputs, and stage commands are bound in the retained
`config.json` and `audit.json`. The experiment used one cold, ordered 3+2
five-summand search call per arm with a 5,000-conflict cap.

| Arm | K600 pool | Shared signed orbits | Residual signed orbits | Search wall L0 | Conflicts | Models | Disposition |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Gray | `public_x_gray_prefix` | 503 | 97 | 265.547 s | 5,000 | 0 | `UNKNOWN_search_cap` |
| Hash | `public_x_hash` | 0 | 600 | 221.433 s | 5,000 | 0 | `UNKNOWN_search_cap` |

Both arms imported the same sequential K600 base. Both Docker builds,
imports, and searches exited cleanly without OOM or an outer timeout. Each
search attempted one identity pattern and one solver call. Neither produced a
candidate. The 5,000-conflict stop is a censored observation, not a solver
refutation. The different residual-domain sizes and one-run walls do not
rank the bases by relation throughput or total index-calculus cost.

The 54-base regeneration took 213.151 seconds of hosted preparation wall;
the paired supervisor observed 930.359 seconds. The conservative reserved
charge is 4,800 seconds for the two arms plus 214 seconds for preparation,
or **5,014 seconds**. This charge is separate from observed wall and does not
include future rank, linear algebra, individual-log, verification, or matched
rho phases. `selected_best_total_runtime` remains `null`.

Independent `audit_replay.py` passed on the downloaded 51-file GitHub
artifact. It verified the three compressed SHA-256/BLAKE3 and plain BLAKE3
hashes against the 54-base pilot manifest; stage stream hashes, clean Docker
states, source and workflow commits, exact signed-orbit partitions, search
configuration, and both raw result dispositions. Changing a result status to
`PASS_one_bounded_partial` without a candidate was rejected for each arm.
`audit.json` contains the SHA-256 and byte count of every retained file.

Fresh S3 GETs of the three canonical factor-base objects matched the hosted
archive bytes exactly; see `s3_roundtrip.json`. The full hosted artifact was
also uploaded under the content-addressed S3 URI in `s3_artifact.json`, then
downloaded with the same SHA-256. Replaying the extracted S3 copy reproduced
`audit.json` byte for byte. The S3 copy preserves the raw evidence beyond
the GitHub artifact's retention period.

Reproduce the retained checks with `blake3==1.0.8` available to Python:

```sh
gh run download 38983406418 --dir /private/tmp/n83-hosted-run-38983406418
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-paired-run-38983406418/audit_replay.py \
  --artifact /private/tmp/n83-hosted-run-38983406418/n83-hosted-pair-38983406418-1 \
  --receipt /private/tmp/n83-retained-audit-replay.json
```

Use a fresh receipt path because the verifier refuses to overwrite evidence.
The next decision requires natural relation yield and useful rank across
multiple matched rounds. Complete-runtime selection additionally requires
the remaining index-calculus phases and matched rho on the same public
targets; neither arm can be selected from this first bounded pair.

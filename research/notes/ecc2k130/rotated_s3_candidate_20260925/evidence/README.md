# Recursive-S3 candidate evidence

`receipt.json` records the complete sequential n13 then n19 producer and
independent-verifier run, with each child command, UTC interval, exit code,
stdout/stderr hash and result hash. The full producer `result.json` under
each arm retains all 32 target classifications and every terminal-matching
candidate path; `verify.json` records the independent replay. The empty
stdout/stderr files are retained because an empty child stream is a measured
fact, not an implicit success or UNSAT marker. `progress.json` preserves the
last completed mask ordinal. The archive is in this Git PR, so its durable
location is the repository path rather than a local temporary directory.

`manifest.json` hashes every evidence file other than itself. The original
n13 JSON-key mismatch and its first freeze are retained in `failure_0`.
The corrected full run used the final `FROZEN.json` SHA-256
`1edd91441f0d17e882c8b5905bf15abacfef1daf2cc435aab210a13e1aa818de`.
The absolute `/private/tmp` paths in `receipt.json` identify the measured
host run; reproduce it on another checkout with a fresh output path:

```sh
python3 research/notes/ecc2k130/rotated_s3_candidate_20260925/ci_replay.py
python3 research/notes/ecc2k130/rotated_s3_candidate_20260925/run.py --out /private/tmp/rotated-s3-candidate-replay
```

For an independent replay of a committed producer output without repeating
the producer, use `verify.py --arm n13-m5 --producer
research/notes/ecc2k130/rotated_s3_candidate_20260925/evidence/n13-m5/producer/result.json
--out /private/tmp/rotated-s3-n13-verify.json`, and analogously for n19.
The committed source pins every inherited input hash and rejects drift.

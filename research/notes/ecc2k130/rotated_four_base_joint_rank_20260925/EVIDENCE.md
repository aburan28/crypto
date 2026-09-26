# Raw evidence and replay

The pre-outcome head `878d37dbb88c44a252dce1892c1c63c4b7c89d51` in draft [PR #779](https://github.com/aburan28/crypto/pull/779) froze the protocol, all code, the literal unconditional 256/64 target split, and SHA-256 values before rank or holdout results. Its hash-only `parse` and `replay` CI jobs passed before the measured run. `FROZEN.json` SHA-256 is `28c122c6ecbf08efc846510f65658ef42092480128e44f8994f1631b67896549`. The source archives are the complete independently replayed #767/#769 censuses, pinned to `39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c` and `fe84aef6a2cf7f6f4c950245c9c8e870354fb750666f5482b81f9a997d107140` respectively. The exact code and input hashes are in `FROZEN.json`; the literal input-construction operation counts and hashes are in `inputs/construction_receipt.json`.

Run from a checkout with Python 3.13:

```sh
python3 research/notes/ecc2k130/rotated_four_base_joint_rank_20260925/ci_replay.py
python3 research/notes/ecc2k130/rotated_four_base_joint_rank_20260925/run.py --evidence /tmp/rotated-joint-new-evidence
python3 research/notes/ecc2k130/rotated_four_base_joint_rank_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_four_base_joint_rank_20260925/evidence
```

The first command verifies the frozen source, archives and input split. The second creates **new** evidence without overwriting this first run. The third verifies SHA-256 of the committed compressed raw archive, extracts it safely, and independently replays the selected rows, exact modular rank, point-only archive lookup, group-checked recovered scalars and sealed labels. CI runs the first and third commands. The committed [raw.tar.gz](evidence/raw.tar.gz) is 26,309 bytes with SHA-256 `00ceb79d014bd92b4ce54baa9c1923f1594532b719881af17ae44a9e1aa1fe98`; its internal files and hashes, all exact child argv/UTC intervals, caps, exit codes, stdout/stderr, CPU/wall and high-water memory are in [receipt.json](evidence/receipt.json). It stores training rows, base logs, archive responses, point-only results and the independent verification result. The prior complete `7^6` tuple censuses are referenced by content hash rather than copied into this tiny diagnostic.

The first run occupied 14:11:38–14:11:47 UTC on 2026-09-25. Training/oracle/recovery/independent-replay children took 3.307/2.929/0.061/2.099 seconds wall, all exit zero under their 120/120/120/300-second caps. Each training and oracle child loaded 23,927,716 uncompressed archive bytes. Their own high-water RSS values were 313,278,464 / 311,279,616 bytes; the recovery and verifier reported 23,281,664 / 312,410,112 bytes. The run retained operation counters in each raw JSON. These are archive-stage host costs and do not price discovering an implicit PDP witness or certifying a negative one.

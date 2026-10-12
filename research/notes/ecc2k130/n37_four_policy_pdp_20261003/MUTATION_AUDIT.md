# Replay rejection controls

Every control below used the preserved pre-Clippy expanded `RESULT.json` with SHA-256
`5efc8fd85e974456ee25a84ce47e373a6cea5cee4b73dc3dbfd998aa6ed8f799`.
The changed result files were scratch files under `/private/tmp`; only their
small failure receipts are retained here. Regenerate a control by expanding
`RESULT.json.gz`, applying the stated `jq -c` expression to the expanded
file, and passing that changed file to `n37_four_policy_pdp_replay` with a
new receipt path. The replay must exit nonzero and write `"status":"FAIL"`.

| Control | Exact scratch transformation | Retained receipt | Expected error |
| --- | --- | --- | --- |
| Unpaired witness | `.policies[0].targets[0].m3.indices[0] += 1` | `MUTATION_WITNESS_REPLAY.json` | paired target rows differ |
| Paired witness | `.policies[0].targets[0].m3.indices[0] += 1 \| .policies[1].targets[0].m3.indices[0] += 1` | `MUTATION_PAIRED_WITNESS_REPLAY.json` | exact target-0 PDP result differs |
| Paired rank row | `.policies[0].rank_trace[0].row[0] += 1 \| .policies[1].rank_trace[0].row[0] += 1` | `MUTATION_PAIRED_ROW_REPLAY.json` | rank probe 0 differs |
| Paired hit flag | `.policies[0].targets[0].m2.status = "hit" \| .policies[1].targets[0].m2.status = "hit"` | `MUTATION_PAIRED_HIT_REPLAY.json` | exact target-0 PDP result differs |
| Decision | `.selection_comparison.decision = "FIXED_BLOCK_M3_SELECTION_LEAD"` | `MUTATION_DECISION_REPLAY.json` | replayed source-native decision differs |

For the public-input control, the first line of
`n37_L1024_b03.points.jsonl` was temporarily changed from
`[46335234820,116989936517]` to `[46335234821,116989936517]` while the
unchanged final manifest was replayed. `MUTATION_TARGET_REPLAY.json` records
failure at the pinned input hash, with mutated SHA-256
`8039962c3fd5f752709339848b53604101c3b66cd876d8babdc412dabf2a8bf6`.
The original file was restored in a shell `EXIT` trap and verified at its
frozen SHA-256
`84400a2914f06e4a001d0f113f0195952f692634d4ff8285e23ade2599e2bde2`.

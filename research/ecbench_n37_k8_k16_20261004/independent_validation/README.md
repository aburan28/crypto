# Independent full replay

[CI run 37191831916](https://github.com/aburan28/crypto/actions/runs/37191831916)
built `ecbench` on Ubuntu 24.04 x86-64 and ran:

```sh
./target/release/ecbench verify \
  --dir research/ecbench_n37_k8_k16_20261004/sessions/mac_arm64_l0_01 \
  --replay-all --exit-code --out RECEIPT.json
```

The [receipt](RECEIPT.json), SHA-256
`8a5f8683c8b3aea4c1e8343242d4d6319f4938fb8f6ae4a342c9c4ad3dad5baa`,
reports 384/384 verified records and 320/320 identical measured replays,
with zero problems. Its auditor environment class
`ECBENV2hda64b23436f4` differs from the Mac producer's
`ECBENV2hd82681268e96`. It verifies deterministic outcomes and counted
work, not an isolated wall-time speedup.

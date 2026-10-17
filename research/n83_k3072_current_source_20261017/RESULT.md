# N83 K3072/root11 works on the current IC source

The current finalist source, SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`,
completed one K3072/root11 fresh process on the previously exposed N83
single-point fixture. Independent arithmetic replay checked all 509,952
distinct factor-base points, the four-point relation, the recovered scalar,
full rank 3,072, and the five exclusive target-online phases. The exact
current-source candidate is
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a`.

| Source commit | Base hash | Root entries | Local cold process | Target online | Replay |
| --- | --- | ---: | ---: | ---: | --- |
| `416489e09c3d8885b6c12f262935e6ba29373460` | `c5511247735bc183359d5425c4dd794aaad49ad38d2ca964799e1caa76afaa08` | 783,270,873 | 678.111 s | 31.908 ms | `PASS_INDEPENDENT_ARITHMETIC_REPLAY` |

The Rust release binary SHA-256 was
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
The producer ran under a 900-second child cap on an Apple M4 Pro with 48 GB
RAM. Local process scheduling was shared; these wall times are feasibility
diagnostics. The base hash and root-entry count match the earlier verified
K3072/root11 Tier 1 candidate, while the exact IC1 identity changes with
the executed source commit and content hashes. The [candidate manifest](candidate_IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a.json),
[producer receipt](receipt.json), and [replay summary](replay.json) retain the
identity and stage outputs.

The complete raw run, including the 40 MB base dump, target record, build
log, and exact local replay scripts, is archived at
`/Volumes/SSD990/crypto/artifacts/n83_k3072_current_source_20261017/raw.tar.gz`,
SHA-256 `ba747b054d24af88abfb6df66acbd2efb407b5c92fb96092b6f2f0b4aa7f65d1`.
An earlier archive with an inherited source-commit field is retained beside
it as `pre_identity_correction.tar.gz`, SHA-256
`e472296e93d626c30e753c3b73d17ae57b03baf259df8f320bd12fa6d8850048`;
the current manifest uses the actual executed repository commit above.
After extracting the archive, the versioned `verify_archive.py` checks its
SHA-256, the source and input hashes, the full base, relation, scalar, phase
accounting, and candidate manifest. It passed against the archived files:

```sh
tar -xzf /Volumes/SSD990/crypto/artifacts/n83_k3072_current_source_20261017/raw.tar.gz -C /path/to/extracted
python3 research/n83_k3072_current_source_20261017/verify_archive.py /Volumes/SSD990/crypto/artifacts/n83_k3072_current_source_20261017/raw.tar.gz /path/to/extracted
```

The next comparison should pair this current-source K3072 candidate with
the winner of the frozen K1104–K1120 finalist panel on fresh public points,
with each IC arm paired to a same-point strong rho run and per-pair Tier 1
calibration. The present point was exposed before that panel and has no role
in its selection.

Ten new public points for that followup are sealed before the finalist result
is inspected. `generate.py --verify-sealed` replays the earlier exclusion sets
and all ten scalar/subgroup certificates. Their public-point SHA-256 is
`9bf9bfca72d4e82a5780759d52ac000b812ca25e831af04238983c8b72d2d777`;
the answer sidecar SHA-256 is
`3f7a75dc0134f883a77fa2b87eb6839b6df804c53950d8dc77e7f29b2500284c`.
An independent read-only inventory scanned 490 `research/n83*/**/*.jsonl`
files containing 32,366 earlier two-string point rows and found no overlap
with these ten points.
The answer sidecar is for independent replay, not a producer input. A run
protocol should choose its second candidate by the preregistered finalist
online selection output. If that selection is null or the panel is
unadmitted, the fixed fallback is K1120/root11
(`IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234`),
chosen from the earlier exposed-point screen before the ten points above
are used. The fallback does not turn a null finalist result into a selected
winner. Use two dedicated Tier 1 host strata of five same-point IC/rho pairs.
The versioned runner, independent verifier, two-metric paired selector,
provider-archive binder, IC1 row materializer and cloud workflow implement
that preregistered next step. The workflow waits for a later `dispatch.json`
bound to the independently replayed finalist result; the ten-point fixture
has not been dispatched.

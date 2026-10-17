# K2816/root11 current-source N83 control

The current finalist IC source built the K2816/root11 dynamic index over
**467,456** distinct subgroup-usable factor-base points and **2,816** folded
columns. It reached full rank, built **658,165,181** exact root entries,
and recovered the exposed public point's scalar with a four-point relation.
Independent arithmetic replay reconstructed every factor-base point and
checked the relation, scalar, five online phases, and exact candidate identity.
The candidate is
`IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b`.
The [local IC1 run row](../n83_smallk_current_source_20261017/local_ic1_rows.json)
uses workload ID `eebf98ee1031` and run ID
`IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8bWeebf98ee1031R1`.

| Current-source exposed-point control | K2543/root11 | K2816/root11 | K3072/root11 |
| --- | ---: | ---: | ---: |
| Actual usable base points | 422,138 | 467,456 | 509,952 |
| Effective columns | 2,543 | 2,816 | 3,072 |
| Exact root entries | 536,737,856 | 658,165,181 | 783,270,873 |
| Packed root allocation, bytes | 5,368,709,120 | 5,368,709,120 | 5,368,709,120 |
| Local receipt wall, s | 401.115 | 431.024 | 679.051 |
| One-target online, ms | 18.983 | 25.301 | 31.908 |
| Independent arithmetic replay | Pass | Pass | Pass |

K2816's setup was **390.719 s**, including a **389.683 s** exact-index
build; its rank stage was **39.796 s**. The one-target interval comprised
0.000 ms query, 25.183458 ms point decomposition, 0.012 ms relation check,
0.000625 ms descent, and 0.105083 ms recovery check. Its target used 79,962
probes. The process cap was 1,000 s and the measured process finished in
431.024 s. The source SHA-256 was
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`;
the local release binary SHA-256 was
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
The base hash was
`9208460bf63be634f0aee6091abaaed855e05004fc2894e28a9c1fd394c35ec4`.

K2816 was 1.333 times K2543 in this exposed-point online interval and
1.075 times its cold process duration. These are shared-host diagnostics on
one exposed point. The earlier K2543 and K3072 values are in their versioned
current-source result directories. The next candidate-selection step uses
unused public points in alternating Tier 1 same-point IC–rho pairs.

Raw output including the 36 MB base dump, command, flags, receipt and local
replay scripts is committed as [raw.tar.gz](raw.tar.gz),
SHA-256
`61daed9917be546f46bf72cc8d634525fefbce5a27a01cf72fd56fc5f6cfde63`.
With that archive extracted to a scratch directory, the
[versioned verifier](verify_archive.py) checks the archive digest and replays
the source, flags, entire base, relation, scalar, phase accounting, and
candidate ID:

```sh
mkdir -p /private/tmp/n83-k2816-replay
tar -xzf research/n83_k2816_current_source_20261017/raw.tar.gz -C /private/tmp/n83-k2816-replay
python3 research/n83_k2816_current_source_20261017/verify_archive.py research/n83_k2816_current_source_20261017/raw.tar.gz /private/tmp/n83-k2816-replay
```

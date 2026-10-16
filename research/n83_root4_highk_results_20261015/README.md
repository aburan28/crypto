# N83 root4 high-K balanced-order records

This directory seals the four-arm K3550/K3596 root4 comparison in
[`PROTOCOL.md`](../n83_root4_highk_20261015/PROTOCOL.md). Candidate manifests
name the exact base and root index. Four measurement rows preserve both fresh
processes per candidate and the common 4,000-point public panel.

The executed source commit is `8a37349f3d3cc2bc2960b03167adcf57da83ba77`.
The public-point SHA-256 is
`48a243619e270abc096e172fc669bb03494fd8c09d4e39c55935121408d44eea`;
the known-scalar sidecar SHA-256 is
`fa660a434e2618a6fb3a9c251843c6bc4c9cc596c612c244ca2bb2ff8d54262d`.
The sidecar is used for independent replay outside the measured intervals.

The sealer checks the exact raw ZIP and independent analysis digests, source,
fixture, base, resource envelope, and analysis identities. Repeating the
command verifies byte equality of each canonical record. The
[result report](RESULT.md) retains the bounded ranking and the checked
[capacity proof](Capacity.lean). Artifact `11967981375` from run `38719753677`
has ZIP SHA-256
`1927c1dcbbdc2ccabaa708896624e5e8735359462fc2edfec12464442cbc158b`;
the independently replayed `analysis.json` SHA-256 is
`a8f7ec556a66212400b629ff9ce85980735b609f0802b9c6e8e3cf2ec811152e`.

```sh
python3 research/n83_k3250_sparse_fresh_20261015/generate.py
python3 research/n83_k3350_sparse_fresh_20261015/generate.py
python3 research/n83_k_abba_20261015/generate.py
python3 research/n83_root4_highk_20261015/generate.py
python3 research/n83_root4_highk_20261015/verify.py /path/to/results
python3 research/n83_root4_highk_results_20261015/materialize.py \
  --artifact-zip /path/to/raw.zip \
  --artifact-sha256 1927c1dcbbdc2ccabaa708896624e5e8735359462fc2edfec12464442cbc158b \
  --analysis-sha256 a8f7ec556a66212400b629ff9ce85980735b609f0802b9c6e8e3cf2ec811152e \
  /path/to/results
(cd research/n83_root4_highk_results_20261015 && lean Capacity.lean)
```

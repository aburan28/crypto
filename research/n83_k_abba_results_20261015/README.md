# N83 balanced-order factor-base records

This directory seals the four-arm comparison in
[`PROTOCOL.md`](../n83_k_abba_20261015/PROTOCOL.md). Each of the two
candidate manifests names an exact factor base and kernel. Four measurement
rows retain both fresh processes per candidate. The workload record fixes
the common 4,000-point public panel and the experiment's process order.

The executed source commit is `d5e6fc717bac432efba63564ba12081d464eebb6`.
Its public-point SHA-256 is
`35e8b43966d15065cf0186ace762f823a602aa1be7d7c0a9039e020550d01491`;
the known-scalar sidecar SHA-256 is
`6e1181ba958708803316903a645463425bc79130074041d39c8094778e2f4379`.
The sidecar enters independent replay after production, outside measured
target intervals.

The sealer requires the exact raw ZIP and independently generated analysis
digests. It checks source, fixture, base, resource, and analysis identities
before writing canonical candidate, workload, measurement, selection, and
index records. Repeating the command verifies byte equality of existing
records. The [result report](RESULT.md) records the bounded comparison and raw
receipts. Download artifact `11964534821` from run `38713270151`, then
verify its ZIP SHA-256
`62f5c35620f993203a6842f1ac643c56b05067db881214d5298177ccfff65bdc`.
The independently replayed `analysis.json` SHA-256 is
`65e66d4c5f3052eadc3dd3c23ca2df3d2203c67d10b08d2c7cba205936f59cde`.
The source-pinned CI workflow regenerates the fixtures, replays the raw
artifact, and requires byte-identical current-head records.

```sh
python3 research/n83_k3250_sparse_fresh_20261015/generate.py
python3 research/n83_k3350_sparse_fresh_20261015/generate.py
python3 research/n83_k_abba_20261015/generate.py
python3 research/n83_k_abba_20261015/verify.py /path/to/results
python3 research/n83_k_abba_results_20261015/materialize.py \
  --artifact-zip /path/to/raw.zip \
  --artifact-sha256 62f5c35620f993203a6842f1ac643c56b05067db881214d5298177ccfff65bdc \
  --analysis-sha256 65e66d4c5f3052eadc3dd3c23ca2df3d2203c67d10b08d2c7cba205936f59cde \
  /path/to/results
```

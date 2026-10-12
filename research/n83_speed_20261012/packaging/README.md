# Evidence archive

`package_evidence.py` packages tracked, already verified files and rendered PDFs.
It performs no research computation or experiment. Run after the final artifacts
and MANIFEST.sha256 are staged:

```sh
python3 research/n83_speed_20261012/packaging/package_evidence.py
```

Output is create-new: /workspace/report_output/n83_faster_rho_evidence_bundle.zip,
with a separate packaging receipt there. CRC is checked after writing. The zip
is not added to its own contents or checksum manifest. The bundle includes this
study, previous N83 evidence, archived exact-curve inputs, measured worker sources,
the Lean CI workflow and the comprehensive70-page PDF. Frozen paths remain
original; reproduction uses a new protocol and never rewrites old records.

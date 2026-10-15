# Final N83 tournament admission

`PROTOCOL.md` froze the candidate rule and 4,000 fresh public points before
the four source screen rankings were available. `admit.py` implements that
rule. It requires all four source ZIPs and independently replayed extracted
directories, plus their versioned result packages, before it writes
`admission.json`.

Arrange the archived sources under one directory as:

```text
n83-best-e2e-38768544528.zip
n83-best-e2e-38768544528/results/analysis.json
n83-legacy-challenger-38770004923.zip
n83-legacy-challenger-38770004923/results/analysis.json
n83-root4-bloom-38772262240.zip
n83-root4-bloom-38772262240/results/analysis.json
n83-root11-highk-38773883285.zip
n83-root11-highk-38773883285/results/analysis.json
```

Each extracted directory also retains its `source.json`, `host.json`, and
complete raw runs. From this repository:

```sh
python3 research/n83_final_tournament_20261015/admit.py \
  --artifact-dir /absolute/path/to/archived-sources
```

The script checks the pinned source commits and fixtures, the result
packages' artifact and analysis digests, each eligible repeated arm, and
the candidate manifest. It adds K3300/root11/wide3 and K3596/root4/off as
controls, takes the lowest eligible mean from each screen, and deduplicates
exact base, kernel, and source configurations. The output is deterministic
and refuses to overwrite a different admission record. Arm execution and
the final 99% paired decision follow `PROTOCOL.md`.

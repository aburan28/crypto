# N83 final tournament result package

The frozen admission in `../n83_final_tournament_20261015/admission.json`
selects four complete factor-base pipelines from four source screens. The
GitHub run numbered 38793258972 tests each pipeline in list and reverse-list
order on 4,000 fresh public subgroup points. Its source commit is
`d72ee1c00da72361c81060a15f1f25fc828a9d0f`.

After the run, archive its raw GitHub artifact and verify the ZIP against the
artifact digest. Extract it outside this checkout. Independently replay all
base constructions, point relations, recovered scalars and resource receipts:

```sh
python3 research/n83_final_tournament_20261015/verify.py \
  /absolute/path/to/extracted/results
```

Compare the local `analysis.json` SHA-256 to the cloud analysis SHA-256. Then
run `materialize.py` with the verified ZIP, analysis hashes and GitHub artifact
ID. The script freezes a compressed 4,000-point workload, all eight run rows,
the full frozen comparison, and an index of raw receipt hashes. It refuses a
different value for any previously sealed file.

Only a `SELECTED` tournament result can populate
`../n83_final_tier1_20261015/selection.json`. That later experiment measures
five previously unseen single targets against same-point rho on one cloud
instance under the Tier 1 resource and noise gate.

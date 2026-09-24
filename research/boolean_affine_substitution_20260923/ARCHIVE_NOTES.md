# Frozen artifact custody

The primary producer's final manifest includes
`run_01/__pycache__/analyze.cpython-313.pyc` (51,890 bytes), SHA-256
`56de81602fcfcd7f0156d712202b4a670ac73e6d3066d910602ed5e376d96d39`.
Importing the frozen analyzer during a reporting smoke check created this
incidental interpreter cache before the producer sealed the directory. It was
not a measurement input.

The cache is preserved byte-for-byte because deleting it would break the original
manifest. It is the sole explicitly tracked cache in this study. Evidence tests
hash it as an archive member; they import the readable top-level analyzer and do
not execute this bytecode. Ordinary top-level interpreter caches remain ignored.

The primary source snapshots, protocol, raw observations, receipts and manifest
are unchanged. The top-level summary, ledger, correctness note and conclusion are
post-run reporting artifacts. The ledger binds the summarizer and the recorded
confirmation refusal. That refusal has an observed exit status and assertion;
its exact execution timestamp was not separately recorded and remains null.

The primary run is rejected. `confirmation_plan.json` records an unexecuted plan;
`confirmation_protocol.json` and `run_02` do not exist. No confirmation timing or
external-review conclusion is inferred from the plan or producer tests.

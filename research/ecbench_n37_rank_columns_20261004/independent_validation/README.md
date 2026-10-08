# Independent all-run replay

[GitHub Actions run 37182953423](https://github.com/aburan28/crypto/actions/runs/37182953423)
checked commit `fc462a299f0479132d4647c50244e65c5b5eea` on Linux x86-64.
Its `ecbench-n37-rank-columns-independent-receipt` artifact is copied
unchanged as [RECEIPT.json](RECEIPT.json), SHA-256
`efb14a5cb1c60291accc4d28518178f6a7d7ed5305178459936252b346b94029`.
The auditor environment class `ECBENV2hda64b23436f4` differs from the
producer's `ECBENV2hd82681268e96`. The receipt verifies the five session
file hashes, all 384 execution records and exact reproduction of all 320
measured outcomes and deterministic counts. This establishes independent
replay of the arithmetic and accounting; it does not promote the producer's
L0 wall timings or price omitted native work.

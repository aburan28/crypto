# Independent full replay

[CI run 37196179536](https://github.com/aburan28/crypto/actions/runs/37196179536)
built the frozen PR revision `170becb5bc485507cd9b2accce3863e85ba8c3d1`
on Ubuntu 24.04 x86-64, then ran `ecbench verify --replay-all --exit-code`
against the committed Mac session. The saved [receipt](RECEIPT.json) reports
384/384 verified native records and exact reproduction of all 320 measured
deterministic executions, with zero problems. The auditor environment class
was `ECBENV2h7797d0021f24`, distinct from the measuring Mac class
`ECBENV2hd82681268e96`. The receipt's SHA-256 is
`b75740f1a28dd4a5352c15e0e738a82bf9c8d3402e3f6ea797d129c424e05472`.

That same job profiled the 64 frozen first-round children under Valgrind
3.22.0 and checked every recovered scalar. The full raw artifact, including
this receipt, is in [RAW-CALLGRIND.zip](../RAW-CALLGRIND.zip). The
[result](../RESULT.md) keeps the Mac L0 timing and Linux instruction count
in their separate units.

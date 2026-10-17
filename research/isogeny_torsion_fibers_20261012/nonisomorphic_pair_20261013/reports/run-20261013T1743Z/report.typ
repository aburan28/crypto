#set page(paper: "us-letter", margin: 0.62in)
#set text(size: 9pt)
= Non-isomorphic pair V1: timing wrapper failure
2026-10-13 17:43 UTC · shared Darwin arm64 · L0

The frozen degree-three map, rational kernel and non-isomorphism preflight
passed. The native child wrote 24 verified rows; the outer `/usr/bin/time -l`
returned exit 1 on a restricted-host `sysctl kern.clockrate` denial.
No wall-time result is admitted. The counted group-operation graph below is
derived from the preserved `raw-v1.tsv`. Full commands, hashes, inputs and
stderr are in `README.md` and the parent directory.

#image("group_ops-v1.svg", width: 100%)

V1 totals across eight targets: negation rho 664, source relations 3,312,
destination relations 3,794 counted group additions and doublings.
The clean direct v2 execution is a separate run and has its own receipt.

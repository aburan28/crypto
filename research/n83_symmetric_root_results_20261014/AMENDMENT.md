# Symmetric-index witness and density observation

The frozen [protocol](../n83_symmetric_root_20261013/PROTOCOL.md) expected
the quotient to preserve the first full-state witness for every canonical
root key. The six-representative N83 source test passed that property on its
fixture. At K1798, the measured symmetric and full-state arms have the same
base hash `e8e34092664026af82eb730ec5c0ccc3e4f249e26d49aae12783f1da8d68cd75`
and the same **count** of 268,316,606 distinct root keys. Their target rows
nevertheless differ: the symmetric arm used 119,686 probes and point indices
`[157782,273832,202015,217163]`; the full arm used 111,004 probes and
`[31597,175234,250462,263363]`. Both four-point relations and the recovered
scalar passed independent replay.

The target-row difference does not isolate a lookup-witness change: rank
worker completion can also change the later query sequence. First-witness
parity at K1798 therefore remains unverified by this receipt. Key-count
parity alone cannot establish it. The large run did not export an ordered
digest of all root keys, so it also does not by itself establish key-set
equality. The Lean
S3/Frobenius lemmas cover the algebraic transport; the source canonicalization
and insertion-order bridge stays open in the proof map. Future descriptions
of the quotient use the replayed target witness and measured key count
without equating the full and quotient first-witness maps.

The same K2543 receipt records 536,737,856 root keys from 268,375,505
retained states. This overturned the initial `ceil(2S/3)` dense-table request
in the separate follow-up branch: it would leave too few slots at larger K.
That branch was revised before push to request `ceil(4S/3)`, preserve the
five-byte exact lookup, and test K2830/K2800/K2750 with a 14-GiB RSS guard.

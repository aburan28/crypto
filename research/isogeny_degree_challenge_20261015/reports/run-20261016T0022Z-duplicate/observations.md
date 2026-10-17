# Duplicate-kernel completeness regression

The native rejection test failed against frozen source
`ff71f577203179760211cdc1dc3da7a33fc77ec1`: the judge returned ACCEPTED, degree
13, coverage COMPLETE after one valid P192/13 isogeny was copied into both map
slots and the construction receipt hash was recomputed. The mathematical map
remained valid. The two slots supplied one distinct kernel, so the completeness
claim was unsupported.

The regression source is commit `2934a2eca`; the corrected runtime is
`9bb7cbce2b84151eff6846b02da842a8ef019343`. The correction compares kernels
modulo the source characteristic and requires distinct Frobenius roots before
accepting COMPLETE coverage. Partial valid maps remain eligible; repeated-root
eigenline totals are left unestablished. The corrected full native package
suite passed, and its installed replay packet is recorded separately.

The run identifier comes from the original command file's creation time in UTC.
The source is the registered P192 model; both records point to the same target.
Exact source/target ICV1 identities, scalar transports and the executed binary
hash are in the [fresh certificate](observed-certificate.json) and
[observed entry](observed-entry.json). The [raw packet](raw.tar.gz) preserves
all mutated bytes, receipt hashes, command outputs and the failing assertion.
The [visual report](report.pdf) embeds the [diagnostic diagram](../../duplicate-flow.svg).

This is a coverage-metadata correctness defect on a public control. It adds no
new discovery or speed measurement and does not change the degree-13 proof.

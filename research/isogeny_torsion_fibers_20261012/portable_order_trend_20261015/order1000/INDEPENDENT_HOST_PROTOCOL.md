# Order-1,201 independent Linux replay

**Frozen before Linux execution.** This experiment runs the selected
`p=1237, r=1201` cell on an independent GitHub `ubuntu-24.04` host. The
source, public and checker inputs, seven arms, target prefixes
`1/4/16/64`, three rotated rounds, relation and rho caps, and 300-second
child limit are fixed by [`RUN_SPEC_V2.json`](RUN_SPEC_V2.json). The
workflow and shell script hashes are fixed by
[`HOST_SPEC_V1.json`](HOST_SPEC_V1.json); changing either requires a
new version and a reason. The independent Rust checker must replay all
84 child records, 1,785 scalar outputs, and 12,486 relation witnesses.

## Question, reference, and boundaries

The question is whether the larger-order counted result survives a
different architecture/OS build and whether its complete CPU cost remains
above the matched rho and BSGS arms on that host. The reference is the
same-target `rho-negation` arm, with `rho-negation-aa` as an A/A noise
control. The measured unit here is process CPU microseconds per cold child.
The already frozen group-operation price uses
`S_group = total counted group operations / sqrt(r)` per target, including
setup, relation work, solve, transport, and scalar verification. This
Linux replay checks the counted vectors exactly; CPU does not replace
`S_group` or provide a conversion to it.

The primary timing summaries are three paired round ratios at prefix
64: full-image/rho, full-image/BSGS, and rho/rho-A/A. Preserve every
round, the median, minimum and maximum, failed child, stderr, exit,
peak RSS, host description, compiler version, source and binary hash.
Failed or timed-out children are censored; an incomplete check does not
enter the ratio summary. The Linux artifact must be copied and SHA-256
indexed in a versioned study directory before reporting a replay result.

Success requires all 84 child exits, the independent checker receipt,
identical per-child counted arithmetic and scalar outputs to the frozen
Darwin run, a host/build identity, and the descriptive paired CPU and
A/A ratios. Any arithmetic mismatch is a correctness investigation;
a timing change without a count change is a host/compiler difference
pending further controls. A failed workflow is preserved with logs and
classified by its actual failure. The hosted VM has no reserved physical
core, so its timing is L1 descriptive and cannot satisfy the separately
requested isolated L2 timing gate. Whole-core NUMA-local execution with
an isolation receipt and admitted A/A gate remains the exit condition for
that gate; the current AWS account has zero Standard On-Demand vCPU quota
in the screened regions, recorded in [`ISOLATION_ACCESS_V1.md`](../ISOLATION_ACCESS_V1.md).

The workflow uploads `isogeny-fiber-order1201-linux-v1`, including every
raw child and its replay receipt. After download, run the native
`different_j_order1201_host_performance` analyzer with the artifact
directory, this study directory, a new JSON output and a new SVG output.
That analyzer checks exact counted vectors against the Darwin run before
writing timing ratios. The final report and benchmark mirror must cite
the Linux workflow URL, artifact SHA-256, analyzer outputs, and any
unmet gate.

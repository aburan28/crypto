# One-shot dispatch amendment

The first workflow commit accidentally allowed the `measure` matrix to run
after every matching pull-request push. This was found before the hosted
validation finished. [Run 36762587606](https://github.com/aburan28/crypto/actions/runs/36762587606)
was cancelled; GitHub reports both `validate` and `measure` as cancelled,
with no executed steps in the measure job. No timed child started and no
cell result exists from that run.

The workflow now runs `measure` only for an explicit `workflow_dispatch`.
Pull requests still run `validate`, including independent replay of the
frozen inputs and local smokes. Dispatch the six-cell matrix exactly once
after validation passes; archive all cell artifacts, including failures.
An outcome-only PR commit may run validation again but cannot repeat the
cold panel. This amendment changes scheduling only: frozen Q, source,
K/prefilter, arms, child limits, estimator and decision gates in
`PROTOCOL.md` remain exactly as registered.

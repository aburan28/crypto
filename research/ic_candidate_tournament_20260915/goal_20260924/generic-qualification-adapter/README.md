# Generic workers in the existing qualification tournament

The prepared optimized incumbent remains the first arm and fixture producer.
Additional registry entries can name `adapter: "generic-v1"` and `generic_build`
pointing to an existing controlled `generic_build.py` output. The tournament
retains that build's executable, source archive, dependency manifest and build
receipt, and requires the same Linux amd64 compiler as the prepared producer.
The generic source is not rebuilt with the prepared producer's Cargo flags.

Generic arms are accepted only in `--qualification` mode with an explicit
`--qualification-protocol FILE`; they cannot enter an improvement round yet.
The existing scheduler, frozen targets, resource limits, process clocks,
uncertainty calculations and replay driver are shared by both adapters. This
avoids a separate comparison engine with subtly different accounting rules.

Before any measured solve, inventory declares the method and arithmetic kernel,
and independently checks the actual usable factor base and folded matrix columns.
It performs no target query or solve. Each measured native/profile report must
match that declaration and pass build, dispatch, query-law, base, matrix, descent,
scalar and exclusive-phase checks. The ordinary and compressed Callgrind parsers
use the same checksum and phase-coverage gates. All generic native/profile
execution IDs are distinct; timeout, OOM and incomplete records retain their
identities and process states, while complete costs and online results stay null.

Primary timing is native one-target online wall time after reusable preparation
through scalar replay. Whole-process instructions and native time are separate
supplementary costs. The paired records do not treat profiled elapsed time as
native timing. Rho has a reference identity and no fictitious IC stages. The
canonical workload identity is common across generic and prepared arms.

The [registered integration control](PROTOCOL.md) exercises the connection on one
toy cell, including intentional collection exhaustion. It is not the five-cell
comparative qualification or an instrumentation-overhead study. The generated
development table remains a control diagnostic; its results cannot alter the
accepted reference binding or authorize promotion. The two remaining bounded
improvement rounds still need calibrated generic/reference comparisons and an
explicit observer-cost treatment before generic candidates may participate.

To run the control on its declared Linux environment, supply a prepared `scaled`
source and a controlled generic build:

```sh
python3.12 research/ic_candidate_tournament_20260915/generic_driver_smoke.py \
  --prepared /tmp/ic-producer --generic-build /tmp/ic-generic-build \
  --out /tmp/ic-mixed-control
```

The output directory must be new. The runner executes the existing frozen
`tournament.py prepare`, `run`, and `verify` sequence and retains raw evidence.
Later confirmation must exclude this control's generated points as well as all
other development exposures; use the existing supplemental-exclusion mechanism.

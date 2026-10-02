# Preserved preparation attempts

## Hosted compile attempt 1

GitHub Actions run `37050163947` used implementation commit
`533745cb8ff50d60d19348aea899b10eaa643302` and did not reach a CUDA
resource conclusion.

- The native job stopped because GCC promotes the existing CUDA-only
  `#pragma unroll` warnings from `packed131.h` under `-Werror`.  The native
  helper itself had already passed locally under Clang.  The additive repair
  retains strict warnings while applying the repository's established
  `-Wno-unknown-pragmas` host flag.
- The CUDA 13.3 container confirmed `nvcc` 13.3.73, then stopped before its
  first build because Git discovery did not cross the workspace mount boundary.
  The additive repair passes the exact checked-out head SHA from the workflow
  and still cross-checks it against Git when Git metadata is discoverable.

Neither failure compiled a candidate kernel, ran a GPU, produced a resource
receipt, or measured throughput.  They are producer-infrastructure failures,
not evidence for or against shared scratch.

## Hosted compile attempt 2

GitHub Actions run `37051703235` used preparation commit
`0c0178ad6538061760707eb225ac89319112d56e`.  All four production clients and
all three device helpers compiled successfully with CUDA 13.3.73 for native
`sm_120`.  The job then stopped in the independent resource auditor because it
treated cuobjdump's function `SHARED` field as function static alone.

The retained artifact establishes the actual distinction:

| arm | registers | ptxas function static | cuobjdump function extent | difference | stack/local/spills |
|---|---:|---:|---:|---:|---|
| control | 126 | 1,792 | 2,816 | 1,024 | zero |
| cache2 | 126 | 22,272 | 23,296 | 1,024 | zero |
| cache3 | 126 | 32,512 | 33,536 | 1,024 | zero |
| cache4 | 126 | 42,752 | 43,776 | 1,024 | zero |

The 1,024-byte difference is the compiled per-block driver reserve already
registered by the static proposal.  The additive audit repair records ptxas
function static, compiled reserve and compiled extent separately.  This
attempt is positive compile evidence but remains superseded because the audit
process exited nonzero and the subsequent source adds the explicit
square-table incompatibility guard.  No GPU or timing operation ran.

## Exact-head compile gate

GitHub Actions run `37053839217`, job `110993636426`, passed at exact source
`fc3941078db32d1a81f2732fd3622c13866208e6`.  The independent audit rebuilt
locally and reproduced `result.json` byte for byte.  Its result SHA-256 is
`6a304c442fc749864acf4cae9fce26c63b61829bd6a1c1788d5094b13d2ae6de`.
The retained compile receipt and resource logs are in
[`results/compile-fc394107/`](results/compile-fc394107/).

This gate compiled native `sm_120` code without a GPU.  The two-block result is
a static capacity calculation from registers and compiled shared extent;
actual driver-reserved bytes, occupancy and helper correctness must still be
read from the target RTX PRO 6000 before timing.

## GPU attempt 1: native preflight failure

The one allocation started from clean merged source
`c942ee0dd547fdaaf4d54ca925741581d0d0ffa6` on the correct RTX PRO 6000 and
CUDA 13.3.73, then exited before compiling or launching a CUDA arm.  The GPU
producer's direct `g++` command for `test_native.cpp` retained `-Werror` but
omitted the repository's `-Wno-unknown-pragmas` compatibility flag.  GCC
therefore rejected existing CUDA-only `#pragma unroll` lines in `packed131.h`.

- Modal app: `ap-MAQ3Rg69gN6UmAqCFApZHs`
- function call: `fc-01M3Z2SAYFR549J7ZM5TXFVV3X`
- recovery token: `fb552f8989c94e3898655ee8f7286fb6`
- GPU UUID: `GPU-08af0821-3ec3-8ce1-36fb-00c4e37e95b5`
- archive SHA-256: `a4fbe337a7e6f80fee271427be5d084da97ccd3389d86d5b4ffc7563803c318c`
- exit code: `1`

This is a producer-infrastructure failure.  It generated no CUDA candidate,
device-helper result, replay, checkpoint, corpus or throughput sample.  The
additive repair changes only that native compile command and adds a source-audit
assertion for the compatibility flag.  A replacement allocation requires root
coordination under the frozen one-panel authorization.

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

# GPU access notes for the isogeny-rho thread

Date: 2026-09-30.

## What was tried

- Local environment has neither `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET` nor
  `RUNPOD_API_KEY`. Packages `modal` and `runpod` are not installed.
- Unauthenticated HEAD probes and the probe script receipt live in
  `access-probe.json`.
- A prior repository receipt
  (`research/notes/ecc2k130/tau_adic_arithmetic_20260929/gpu/results/runpod-access-01.json`)
  already recorded authenticated RunPod inventory as HTTP 403.

## What a GPU run could and could not show

The `ecc2k130/` CUDA campaign is generated for the Koblitz challenge model
only (`generated/eccF131.h` and the ONB/σ walk). It does not accept an
arbitrary isogenous `a₆`. Building a leaf GPU rho is a separate engineering
project.

Even with credentials, a Koblitz-only throughput smoke is a practicality
note for `E_0`. It cannot establish that an isogenous leaf is faster: the
eligible automorphism order on descending leaves is `A = 2`, so expected
iterations are `√131 ≈ 11.45×` larger than on `E_0` before any per-step
constant.

## CI path

Repository Actions already carry Modal secrets for
`.github/workflows/nist-rtxpro6000-run.yml`. This thread does **not**
dispatch a new paid Modal job: the CPU automorphism screen and the derived
`n = 131` floor answer the research question without burning GPU budget on
an unimplemented leaf kernel.

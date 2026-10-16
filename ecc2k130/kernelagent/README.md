# KernelAgent on the ECC2K-130 walk step

This directory is a harness that points [KernelAgent](https://github.com/aburan28/KernelAgent)
(an LLM agent that writes and checks Triton GPU kernels) at one step of the
ECC2K-130 rho walk. The model behind it is hosted on Fireworks AI (GLM 5.3
Flash by default). The question is whether the agent finds an implementation
of the walk step worth porting into the production client.

| File | What it does |
|---|---|
| `walk_step_problem.py` | The target, in KernelBench format: `P ← P + σʲ(P)`, `j = 3 + ((HW(x) >> 1) & 7)`, over GF(2^131) in the permuted type-II normal basis, with the device's special-case-free affine addition. Batched `(N, 3)` int64 limbs in, `(N, 2, 3)` out. |
| `check_reference.py` + `.cpp` | Checks the target bit for bit against `include/ref.h` (`Ref<CfgF131>`, `addPtRaw`). Run it after any edit to the target. |
| `run_kernelagent.py` | Adds Fireworks as a KernelAgent provider at run time (the fork is unmodified), routes every stage to one model and hands the agent a strict exact-match test in place of its float-tolerance default. |
| `bench_candidate.py` | Re-checks the generated kernel on a fresh seed at three batch sizes, then times it against the reference, interleaved, with CUDA events. |
| `modal_kernelagent.py` | Runs both scripts on a Modal GPU (default RTX PRO 6000). Results go to the `ecc2k130-jobs` Volume and come back to `--out`. |

## Run

```bash
pip install modal torch                       # local client; torch only for check_reference.py
python3 check_reference.py                    # target == include/ref.h
export FIREWORKS_API_KEY=...                  # never committed; passed to Modal as a secret
modal run --detach modal_kernelagent.py --out /tmp/ka-walk [--model <fireworks id>] [--gpu B200]
modal run modal_kernelagent.py --fetch-token <token> --out /tmp/ka-walk   # FIREWORKS_API_KEY must still be set (any value)
```

Without `--model`, the run picks the first Fireworks model id that matches
`FIREWORKS_MODEL_MATCH` (default `glm.*(5p3|5.3|5-3|53).*flash`) and records
it in `launch.json`.

## What a result means

`bench.json` reports speedups over a **readable PyTorch reference**, not over
the production client. The production kernel runs a different iteration (the
x-only doubling map with a σ³ bridge, 12.5 B complete scalar updates/s on the
G7; see `../README.md`). A generated kernel is a candidate only. Before any
speed claim it must be ported into the packed CUDA client and measured there
with the paired protocol (`ecbench-measure`, `research-visuals`, and a row in
`benchmarks/README.md`).

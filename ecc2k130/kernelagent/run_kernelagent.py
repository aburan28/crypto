#!/usr/bin/env python3
"""Run KernelAgent on the ECC2K-130 walk step with a Fireworks-hosted model.

Runs inside the GPU container (modal_kernelagent.py), from a KernelAgent
checkout on PYTHONPATH:

    FIREWORKS_API_KEY=... python3 run_kernelagent.py --out /results [--model ID]

KernelAgent ships OpenAI, Anthropic and relay providers only.  This script
registers Fireworks' OpenAI-compatible endpoint as one more provider at run
time, so the fork stays unmodified, and points every stage (router, extract,
dispatch, compose, kernel workers) at one model.  --model takes a Fireworks
model id; without it the id is chosen from Fireworks' live model list
(FIREWORKS_MODEL_MATCH, default GLM 5.3 Flash) and the choice is recorded.

KernelAgent is handed a strict test (TEST_TAIL below): every output limb must
equal walk_step_problem.Model on random and edge-case inputs.  Its usual
LLM-written test would compare with float tolerances, which mean nothing
over GF(2).  The generated kernel lands in <out>/kernel.py, the agent's JSON
summary in <out>/kernelagent.json; bench_candidate.py then times it.
"""
import argparse
import json
import os
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
FIREWORKS_URL = "https://api.fireworks.ai/inference/v1"
DEFAULT_MATCH = r"glm.*(5p3|5\.3|5-3|53).*flash"

TEST_TAIL = '''

# ---- strict test supplied by ecc2k130/kernelagent/run_kernelagent.py -------
def _edge_cases(device):
    rows = [[0, 0, 0], [-1, -1, 7]]
    for b in (0, 1, 63, 64, 127, 128, 130):
        limb = [0, 0, 0]
        limb[b // 64] = 1 << (b % 64)
        if limb[b // 64] >= 1 << 63:
            limb[b // 64] -= 1 << 64
        rows.append(limb)
    return torch.tensor(rows, dtype=torch.long, device=device)


def test_kernel():
    """Every output limb must equal the reference walk step exactly (GF(2) arithmetic: no tolerance)."""
    try:
        from kernel import kernel_function
        device = "cuda"
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA not available")
        ref = Model().to(device)
        g = torch.Generator().manual_seed(12345)
        e = _edge_cases("cpu")
        x = torch.cat([e.repeat_interleave(len(e), 0), _random_elements(BATCH, g)]).to(device)
        y = torch.cat([e.repeat(len(e), 1), _random_elements(BATCH, g)]).to(device)
        expected = ref(x, y)
        result = kernel_function(x, y)
        if not isinstance(result, torch.Tensor) or result.shape != expected.shape or result.dtype != torch.long:
            print("bad output: want int64 %s, got %r" % (tuple(expected.shape),
                  (getattr(result, "dtype", None), tuple(getattr(result, "shape", ())))))
            return False
        bad = (result != expected).any(-1).any(-1)
        if bool(bad.any()):
            i = int(bad.nonzero()[0])
            print("MISMATCH in %d of %d rows; first row %d: x=%s y=%s want=%s got=%s" % (
                int(bad.sum()), len(x), i, x[i].tolist(), y[i].tolist(),
                expected[i].tolist(), result[i].tolist()))
            return False
        return True
    except Exception as exc:
        print("Test failed: %r" % (exc,))
        return False


if __name__ == "__main__":
    sys.exit(0 if test_kernel() else 1)
'''


def register_fireworks(model_hint: str | None) -> str:
    from utils.providers import available_models, models
    from utils.providers.model_config import ModelConfig
    from utils.providers.openai_base import OpenAICompatibleProvider

    class FireworksProvider(OpenAICompatibleProvider):
        def __init__(self):
            super().__init__(api_key_env="FIREWORKS_API_KEY", base_url=FIREWORKS_URL)

        @property
        def name(self) -> str:
            return "fireworks"

        def get_max_tokens_limit(self, model_name: str) -> int:
            return 32000

    model = model_hint
    if not model:
        provider = FireworksProvider()
        if not provider.is_available():
            sys.exit("FIREWORKS_API_KEY is not set")
        ids = sorted(m.id for m in provider.client.models.list())
        pattern = re.compile(os.environ.get("FIREWORKS_MODEL_MATCH", DEFAULT_MATCH), re.I)
        hits = [i for i in ids if pattern.search(i)]
        if not hits:
            glm = [i for i in ids if "glm" in i.lower()]
            sys.exit("no Fireworks model matches %r; GLM models offered: %s"
                     % (pattern.pattern, ", ".join(glm) or "none"))
        model = hits[0]
    available_models.AVAILABLE_MODELS.append(
        ModelConfig(name=model, provider_classes=[FireworksProvider],
                    description="Fireworks-hosted model for the ECC2K-130 walk step"))
    models._model_name_to_config = None
    return model


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default=os.environ.get("FIREWORKS_MODEL"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--rounds", type=int, default=10)
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    model = register_fireworks(args.model)
    print("model:", model, flush=True)

    problem = HERE / "walk_step_problem.py"
    test = out / "test_walk_step.py"
    test.write_text("import sys\n" + problem.read_text() + TEST_TAIL)

    from Fuser import auto_agent
    argv = ["--problem", str(problem), "--test-paths", str(test),
            "--ka-model", model, "--router-model", model, "--extract-model", model,
            "--dispatch-model", model, "--compose-model", model,
            "--ka-workers", str(args.workers), "--ka-rounds", str(args.rounds),
            "--ignore-router-config", "--no-router-cache", "--verify"]
    (out / "launch.json").write_text(json.dumps({"model": model, "argv": argv}, indent=2))

    # auto_agent prints its JSON summary last on stdout; keep a copy of all of it
    import contextlib
    import io
    buf = io.StringIO()

    class Tee(io.TextIOBase):
        def write(self, s):
            sys.__stdout__.write(s)
            buf.write(s)
            return len(s)

        def flush(self):
            sys.__stdout__.flush()

    with contextlib.redirect_stdout(Tee()):
        rc = auto_agent.main(argv)
    text = buf.getvalue()
    start = text.rfind('\n{\n  "route"')
    summary = json.loads(text[start:]) if start >= 0 else {"success": False, "raw_tail": text[-4000:]}
    code = summary.pop("kernel_code", None)
    (out / "kernelagent.json").write_text(json.dumps(summary, indent=2, default=str))
    if code:
        (out / "kernel.py").write_text(code)
    print("kernelagent exit %d, kernel %s" % (rc, "written" if code else "missing"), flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

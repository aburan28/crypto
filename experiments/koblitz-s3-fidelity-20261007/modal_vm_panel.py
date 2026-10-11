"""Run the frozen paired-S3 panel in one Modal VM Sandbox.

The VM exposes more guest diagnostics than the Function container, but guest
CPU pinning, NUMA binding and perf do not prove host-wide exclusive CPUs.
Every output remains exploratory unless a separate strict host receipt exists.

Example:
  python3 experiments/koblitz-s3-fidelity-20261007/modal_vm_panel.py \
      --count 12 --output experiments/koblitz-s3-fidelity-20261007/modal/new_vm_pilot.json.gz

Use a temporary Modal Volume to transfer the full raw ledger. Sandbox stdout
may truncate large payloads. The Volume is removed only after size, SHA-256,
gzip and expected-run-count checks pass; a failed transfer leaves it for
forensic recovery.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import uuid

import modal


HERE = Path(__file__).resolve().parent
APP_SOURCE = HERE / "modal_app.py"
CONSTANTS = {
    "REMOTE_ROOT", "PILOT", "WARMUP", "BIN", "SOURCE_COMMIT", "WORKERS",
    "COLUMNS", "RELATION_SEED", "CPU_REQUEST", "MEMORY_MIB",
}
CPU_OVERRIDE = """
_old_cpu_selection = _cpu_selection
def _cpu_selection():
    result = _old_cpu_selection()
    selected = [cpu for cpu in result['allowed'] if cpu != 0][:WORKERS + 1]
    if len(selected) != WORKERS + 1:
        raise RuntimeError('not enough nonzero guest vCPUs')
    result['selected'] = selected
    result['sampler_cpu'] = next((cpu for cpu in result['allowed'] if cpu not in selected and cpu != 0), None)
    result['affinity_probe'] = _command(['taskset', '-c', ','.join(map(str, selected)), 'true'])
    result['affinity_supported'] = result['affinity_probe'].get('exit_code') == 0
    result['guest_cpu0_excluded'] = True
    return result
"""


def remote_program() -> str:
    """Freeze exactly the stdlib probe/runner functions from modal_app.py."""
    source = ast.parse(APP_SOURCE.read_text())
    kept = []
    for node in source.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            if isinstance(node, ast.Import) and any(alias.name == "modal" for alias in node.names):
                continue
            kept.append(node)
        elif (isinstance(node, ast.Assign) and len(node.targets) == 1
              and isinstance(node.targets[0], ast.Name)
              and node.targets[0].id in CONSTANTS):
            kept.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name.startswith("_"):
            kept.append(node)
    return ast.unparse(ast.Module(body=kept, type_ignores=[])) + CPU_OVERRIDE


def remote_command(program: str, count: int, aa: bool) -> str:
    return (
        "from pathlib import Path\n"
        "__file__ = '/tmp/s3_modal_vm_panel_source.py'\n"
        "Path(__file__).write_text(" + repr(program) + ")\n"
        "exec(compile(Path(__file__).read_text(), __file__, 'exec'))\n"
        "import contextlib, sys, json, gzip, hashlib\n"
        "with contextlib.redirect_stdout(sys.stderr):\n"
        f"    doc = _panel({count}, {aa})\n"
        "doc['kind'] = 'modal_vm_s3_two_curve_panel_v1'\n"
        "doc['vm_guest_topology_is_not_host_isolation_evidence'] = True\n"
        "compressed = gzip.compress(json.dumps(doc, separators=(',', ':')).encode())\n"
        "Path('/results/panel.json.gz').write_bytes(compressed)\n"
        "print(json.dumps({'bytes':len(compressed),'sha256':hashlib.sha256(compressed).hexdigest(),'runs':len(doc['runs']),'blocks':len(doc['blocks']),'summary':doc['summary']}),flush=True)\n"
    )


def run(output: Path, count: int, region: str | None) -> None:
    if output.exists():
        raise FileExistsError(output)
    if count not in (1, 12):
        raise ValueError("count must be 1 (smoke) or 12 (frozen pilot)")
    aa = count == 12
    program = remote_program()
    spec = importlib.util.spec_from_file_location("s3_modal_app", APP_SOURCE)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    app = modal.App("s3-pair-vm-panel-diagnostic")
    volume_name = f"s3-vm-panel-{uuid.uuid4().hex[:20]}"
    volume = modal.Volume.from_name(volume_name, create_if_missing=True)
    print(json.dumps({"volume_name": volume_name,
                      "remote_runner_sha256": hashlib.sha256(program.encode()).hexdigest()}), flush=True)
    with app.run():
        sandbox = modal.Sandbox.create(
            "python3", "-c", remote_command(program, count, aa),
            app=app, image=module.image, cpu=16, memory=8192, timeout=1800,
            region=region, volumes={"/results": volume},
            experimental_options={"vm_runtime": True},
        )
        print(json.dumps({"sandbox_id": sandbox.object_id}), flush=True)
        try:
            sandbox.wait()
            message = sandbox.stdout.read().strip()
            stderr = sandbox.stderr.read()
            if not message:
                raise RuntimeError("VM emitted no panel summary: " + stderr[-6000:])
            expected = json.loads(message)
        finally:
            # A final Volume commit occurs when the Sandbox terminates.
            sandbox.terminate()
        compressed = b"".join(volume.read_file("panel.json.gz"))
        digest = hashlib.sha256(compressed).hexdigest()
        if len(compressed) != expected["bytes"] or digest != expected["sha256"]:
            raise RuntimeError(f"Volume transfer mismatch: {len(compressed)} bytes, sha {digest}")
        panel = json.loads(gzip.decompress(compressed))
        if len(panel["runs"]) != (2 + count * 2 * (6 if aa else 4)):
            raise RuntimeError("incomplete run ledger")
        if len(panel["blocks"]) != 2 * count:
            raise RuntimeError("incomplete target blocks")
        if panel["runner_sha256"] != hashlib.sha256(program.encode()).hexdigest():
            raise RuntimeError("remote runner hash mismatch")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(compressed)
        print(json.dumps({"output": str(output), "sha256": digest,
                          "runs": len(panel["runs"]), "blocks": len(panel["blocks"]),
                          "summary": panel["summary"]}, sort_keys=True), flush=True)
        volume.remove_file("panel.json.gz")
    modal.Volume.delete(volume_name)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--count", type=int, choices=(1, 12), default=12)
    parser.add_argument("--region", default=None)
    args = parser.parse_args()
    run(args.output, args.count, args.region)


if __name__ == "__main__":
    main()

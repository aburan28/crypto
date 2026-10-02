"""A real worker in-process against a real hub: submit, run, measure, collect, read back.

Unprivileged everywhere (direct backend, best-effort policy); the privileged
and podman variants run when the environment allows and prove the cgroup
partition, eviction and container paths.
"""
import asyncio
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

from isolab import protocol
from isolab.blobs import sha256_file
from isolab.fabric import Fabric
from isolab.worker import Worker

PRIVILEGED = os.environ.get("ISOLAB_TEST_PRIVILEGED") == "1" and hasattr(os, "geteuid") and os.geteuid() == 0 \
    and platform.system() == "Linux"


def _podman_ok() -> bool:
    exe = shutil.which("podman")
    if not exe:
        return False
    return subprocess.run([exe, "info", "--format", "{{.Host.OS}}"], capture_output=True, timeout=60).returncode == 0


PY_OK = (
    "import json, os, sys, time\n"
    "out = os.environ['ISOLAB_OUTPUT_DIR']\n"
    "x = sum(i * i for i in range(300000))\n"
    "json.dump({'answer': x, 'repeat': os.environ['ISOLAB_REPEAT'], 'warmup': os.environ['ISOLAB_WARMUP']}, open(os.path.join(out, 'metrics.json'), 'w'))\n"
    "open(os.path.join(out, 'table.csv'), 'w').write('a,b\\n1,2\\n')\n"
    "print('hello from', os.getcwd(), open('params.txt').read().strip())\n"
)


@pytest.fixture
async def lab(nats_url, fabric, tmp_path):
    """A worker on this host, direct backend, in the test's namespace."""
    wf = Fabric(nats_url, namespace=fabric.ns, name="test-worker", lease_s=6.0)
    wf.decline_delay_s = 0.3
    await wf.connect()
    w = Worker(wf, worker_id="w-test", pools=["default"], labels={"site": "test"}, lab_cpus=[],
               state_dir=tmp_path / "state", backend="direct", lock_path=str(tmp_path / "lock"),
               keep_job_dirs=False, block_s=1.0, heartbeat_s=2.0)
    await w.start()
    task = asyncio.create_task(w.run_forever())
    yield w
    w.request_shutdown()
    await asyncio.wait_for(task, 60)
    await wf.close()


def spec(argv, **over):
    s = {"schema": protocol.SPEC_SCHEMA_ID, "runtime": {"backend": "direct"},
         "inputs": [{"path": "params.txt", "content": "seed=7\n"}, {"path": "run.py", "content": PY_OK}],
         "command": {"argv": [sys.executable, *argv]}, "fidelity": {"policy": "best_effort", "settle_s": 0.2},
         "limits": {"timeout_s": 60}, "measure": {"sample_period_s": 0.2}}
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(s.get(k), dict):
            s[k] = {**s[k], **v}
        else:
            s[k] = v
    return s


async def _finish(fabric, job_id, timeout=120):
    rec = await fabric.wait(job_id, timeout)
    assert rec["state"] in protocol.TERMINAL_STATES, f"{job_id} still {rec['state']}: {rec.get('error')}"
    res = await fabric.get_result(job_id)
    assert res is not None, rec
    protocol.validate_result(res)
    return rec, res


async def test_roster_and_submission_validation(fabric, lab):
    ws = await fabric.workers()
    assert [w["id"] for w in ws] == ["w-test"] and "direct" in ws[0]["backends"] and ws[0]["max_tier"] in "ABCD"
    assert ws[0]["default_backend"] == "direct"
    assert ws[0]["calibration"] is None or ws[0]["calibration"]["sha256"]


async def test_success_with_repeats_metrics_and_artifacts(fabric, lab, tmp_path):
    out = await fabric.submit(spec(["run.py"], measure={"warmups": 1, "repeats": 3}, labels={"exp": "e2e"}))
    rec, res = await _finish(fabric, out["job_id"])
    assert res["status"] == "succeeded" and res["outcome_class"] == "completed", res.get("error")
    assert res["worker"]["id"] == "w-test" and res["attempt"] == 1 and res["fence"] == rec["fence"]
    runs = res["runs"]
    assert len(runs) == 4 and runs[0]["warmup"] and not runs[1]["warmup"]
    assert all(r["exit_code"] == 0 and r["wall_s"] > 0 for r in runs)
    assert runs[1]["metrics"]["answer"] == sum(i * i for i in range(300000)) and runs[1]["metrics"]["repeat"] == "1"
    assert res["summary"]["n"] == 3 and res["summary"]["wall_s"]["min"] > 0
    assert res["fidelity"]["grade"] in ("A", "B", "C", "D") and res["fidelity"]["policy"] == "best_effort"
    names = {a["path"] for a in res["artifacts"]}
    assert {"run-1/metrics.json", "run-1/table.csv", "_isolab_logs/stdout.log", "run-1/_isolab/launch.json"} <= names
    assert all(a["stored"] for a in res["artifacts"])
    got = await fabric.get_artifact(out["job_id"], "run-2/table.csv", tmp_path / "t.csv")
    assert got.read_text() == "a,b\n1,2\n"
    entry = next(a for a in res["artifacts"] if a["path"] == "run-2/table.csv")
    assert sha256_file(got) == entry["sha256"]
    log = await fabric.get_artifact_bytes(out["job_id"], "_isolab_logs/stdout.log")
    assert b"hello from" in log and b"seed=7" in log
    assert res["inputs"][0]["kind"] == "content" and res["runtime"]["backend"] == "direct"
    assert res["timing"]["phases"]["measure"] > 0 and "quiesce" in res["timing"]["phases"]
    assert res["placement"]["cpus"] and res["placement"]["isolation"]["lock"]
    checks = {c["name"] for c in res["fidelity"]["checks"]}
    assert "isolation_tier" in checks and "settle_other_cpu" in checks
    assert res["labels"] == {"exp": "e2e"}


async def test_blob_input_round_trip(fabric, lab, tmp_path):
    tool = tmp_path / "tool.py"
    tool.write_text("import os, sys; print('tool ok'); open(os.path.join(os.environ['ISOLAB_OUTPUT_DIR'], 'x'), 'w').write('y')\n")
    digest = sha256_file(tool)
    await fabric.put_blob(digest, tool)
    s = spec(["bin/tool.py"], inputs=[{"path": "bin/tool.py", "sha256": digest, "bytes": tool.stat().st_size, "mode": "0755"}])
    out = await fabric.submit(s)
    _, res = await _finish(fabric, out["job_id"])
    assert res["status"] == "succeeded", res["error"]
    assert res["inputs"][0]["kind"] == "blob" and res["inputs"][0]["sha256"] == digest
    assert "run-0/x" in {a["path"] for a in res["artifacts"]}


async def test_failure_timeout_and_build_step(fabric, lab):
    f = await fabric.submit(spec(["-c", "import sys; print('boom', file=sys.stderr); sys.exit(3)"]))
    _, res = await _finish(fabric, f["job_id"])
    assert res["status"] == "failed" and res["outcome_class"] == "completed" and "exited 3" in res["error"]
    assert res["runs"][0]["exit_code"] == 3 and "boom" in res["logs"]["stderr_tail"]
    t = await fabric.submit(spec(["-c", "import time; time.sleep(30)"], limits={"timeout_s": 1.5}))
    _, res = await _finish(fabric, t["job_id"])
    assert res["status"] == "timeout" and res["outcome_class"] == "not_completed" and res["runs"][0]["timed_out"]
    b = await fabric.submit(spec(["run.py"], build=[{"argv": [sys.executable, "-c", "open('built.txt','w').write('ok')"]},
                                                  {"argv": [sys.executable, "-c", "import sys; sys.exit(9)"]}]))
    _, res = await _finish(fabric, b["job_id"])
    assert res["status"] == "failed" and "build step 1" in res["error"] and len(res["build"]) == 2 and res["runs"] == []


async def test_cancel_running_job(fabric, lab):
    out = await fabric.submit(spec(["-c", "import time; print('started', flush=True); time.sleep(60)"], limits={"timeout_s": 120}))
    for _ in range(100):
        rec = await fabric.get_job(out["job_id"])
        if rec["state"] == "running":
            break
        await asyncio.sleep(0.2)
    assert rec["state"] == "running"
    t0 = time.monotonic()
    assert await fabric.cancel(out["job_id"]) == "running"
    _, res = await _finish(fabric, out["job_id"], 60)
    assert res["status"] == "cancelled" and time.monotonic() - t0 < 30
    assert res["runs"][0]["exit_code"] is None or res["runs"][0]["signal"] or res["runs"][0]["exit_code"] != 0


async def test_unplaceable_job_waits_with_reasons(fabric, lab):
    out = await fabric.submit(spec(["run.py"], resources={"cpus": 4096}))
    for _ in range(60):
        rec = await fabric.get_job(out["job_id"])
        if rec["declines"]:
            break
        await asyncio.sleep(0.25)
    assert rec["state"] == "queued" and rec["declines"][0]["worker"] == "w-test"
    assert "4096" in rec["declines"][0]["reasons"][0]
    assert await fabric.cancel(out["job_id"]) == "cancelled"


async def test_verify_argv_runs_after_timing(fabric, lab):
    s = spec(["-c", "import json, os; json.dump({'k': 17}, open(os.path.join(os.environ['ISOLAB_OUTPUT_DIR'], 'certificate.json'), 'w'))"],
             verify={"argv": [sys.executable, "-c", "import json, os, sys; c = json.load(open(os.environ['ISOLAB_CERTIFICATE'])); sys.exit(0 if c['k'] == 17 else 1)"]},
             measure={"repeats": 2})
    out = await fabric.submit(s)
    _, res = await _finish(fabric, out["job_id"])
    assert res["status"] == "succeeded", res["error"]
    assert res["verification"] == {"status": "verified", "counts": {"verified": 2, "refuted": 0, "no_claim": 0, "error": 0}}
    assert res["runs"][0]["verification"]["status"] == "verified"


async def test_mcp_tools_end_to_end(fabric, lab, tmp_path, monkeypatch):
    import isolab.mcp_server as m
    m._fabric = fabric
    over = await m.isolab_overview()
    assert over["workers"][0]["id"] == "w-test"
    script = tmp_path / "job.py"
    script.write_text(PY_OK)
    sub = await m.isolab_run(command=[sys.executable, "job.py"], backend="direct", policy="best_effort",
                             inputs=[{"path": "job.py", "local_file": str(script)}, {"path": "params.txt", "content": "seed=9\n"}],
                             repeats=2, name="mcp-e2e", labels={"via": "mcp"})
    assert sub["eligible_workers"] == ["w-test"] and sub["uploaded"][0]["sha256"] == sha256_file(script)
    st = await m.isolab_wait(sub["job_id"], 120)
    assert st["state"] == "succeeded", st
    summ = await m.isolab_result(sub["job_id"])
    assert summ["status"] == "succeeded" and summ["summary"]["n"] == 2 and summ["fidelity"]["grade"]
    fid = await m.isolab_result(sub["job_id"], "fidelity")
    assert fid["per_run_checks"][0]["checks"]
    logs = await m.isolab_logs(sub["job_id"])
    assert "seed=9" in logs["stdout"]
    fetched = await m.isolab_fetch(sub["job_id"], "run-1/metrics.json", str(tmp_path / "dl"))
    assert fetched["files"][0]["sha256_ok"] and Path(fetched["files"][0]["path"]).exists()
    jobs = await m.isolab_jobs(labels={"via": "mcp"})
    assert jobs[0]["job_id"] == sub["job_id"]
    with pytest.raises(ValueError, match="no online worker"):
        await m.isolab_run(command=["x"], backend="direct", cpus=4096)
    schema = await m.isolab_schema()
    assert schema["job_spec"]["$id"] == "isolab.job/v1"
    m._fabric = None


async def test_calibration_when_available(fabric, lab):
    import isolab.mcp_server as m
    ws = await fabric.workers()
    if not ws[0].get("calibration"):
        pytest.skip("no C compiler on this host")
    m._fabric = fabric
    sub = await m.isolab_calibrate(worker="w-test", repeats=2, policy="best_effort", iterations=20_000_000, backend="direct")
    st = await m.isolab_wait(sub["job_id"], 120)
    assert st["state"] == "succeeded", st
    res = await fabric.get_result(sub["job_id"])
    checks = {r["metrics"]["checksum"] for r in res["runs"] if not r["warmup"]}
    assert len(checks) == 1  # identical work every repeat
    assert res["summary"]["wall_s"]["cv"] >= 0
    m._fabric = None


@pytest.mark.skipif(not PRIVILEGED, reason="needs root on Linux with ISOLAB_TEST_PRIVILEGED=1")
async def test_privileged_partition_and_eviction(fabric, lab):
    out = await fabric.submit(spec(["run.py"], fidelity={"policy": "standard", "settle_s": 0.5},
                                   measure={"repeats": 2}, resources={"cpus": 1, "smt": "allow", "numa_node": "any"}))
    _, res = await _finish(fabric, out["job_id"])
    assert res["status"] == "succeeded", res["error"]
    iso = res["placement"]["isolation"]
    assert iso["cgroup"] and iso["partition_state"] in ("isolated", "root", "member")
    assert iso["evicted"] is True and iso["irqs_moved"] >= 0
    assert res["fidelity"]["tier"] in ("A", "B")
    assert res["runs"][1]["cgroup"]["usage_s"] > 0 and res["runs"][1]["cgroup"]["memory_peak_bytes"] > 0
    assert res["runs"][1]["timing_source"] == "launcher"
    names = {c["name"]: c for c in res["fidelity"]["checks"]}
    assert names["cpu_partition"]["status"] in ("pass", "info") and names["evicted"]["status"] == "pass"


@pytest.mark.skipif(not _podman_ok(), reason="podman not usable here")
async def test_podman_container_run(nats_url, fabric, tmp_path):
    wf = Fabric(nats_url, namespace=fabric.ns, name="test-worker-podman", lease_s=10.0)
    await wf.connect()
    w = Worker(wf, worker_id="w-podman", pools=["default"], labels={}, lab_cpus=[], state_dir=tmp_path / "state",
               backend="podman", default_image="docker.io/library/python:3.12-slim", lock_path=str(tmp_path / "lock"),
               block_s=1.0, heartbeat_s=2.0)
    await w.start()
    if w.launcher is None:
        await wf.close()
        pytest.skip("no C compiler for the launcher")
    task = asyncio.create_task(w.run_forever())
    try:
        s = spec(["run.py"], runtime={"backend": "podman", "image": "docker.io/library/python:3.12-slim"},
                 measure={"repeats": 2}, placement={"require_images": False})
        s["command"]["argv"] = ["python3", "run.py"]
        out = await fabric.submit(s)
        _, res = await _finish(fabric, out["job_id"], 600)
        assert res["status"] == "succeeded", (res["error"], res.get("logs"))
        assert res["runtime"]["backend"] == "podman" and res["runtime"]["container"]["id"]
        assert res["runs"][0]["timing_source"] == "launcher" and res["runs"][0]["rusage"]["wall_s"] > 0
        assert res["runs"][1]["metrics"]["answer"] == sum(i * i for i in range(300000))
        assert "run-1/table.csv" in {a["path"] for a in res["artifacts"]}
        if shutil.which("runsc") and os.geteuid() == 0:
            s2 = dict(s, runtime={**s["runtime"], "oci_runtime": "runsc"})
            out2 = await fabric.submit(s2)
            _, res2 = await _finish(fabric, out2["job_id"], 600)
            assert res2["status"] == "succeeded", res2["error"]
            assert res2["runtime"]["oci_runtime"] == "runsc"
    finally:
        w.request_shutdown()
        await asyncio.wait_for(task, 120)
        await wf.close()

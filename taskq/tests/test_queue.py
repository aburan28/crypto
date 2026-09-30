import threading
import time

import pytest

from taskq.store import FencedError
from taskq.worker import Worker
from conftest import make_spec


def worker(store, repos, tmp_path, **kw):
    kw.setdefault("lease_seconds", 2)
    return Worker(store, ["cpu"], repos, worker_id=kw.pop("wid", "w1"),
                  artifact_dir=tmp_path / "artifacts", result_dir=tmp_path / "results",
                  block_seconds=0.2, **kw)


def test_command_roundtrip(store, repos, repo, tmp_path):
    tid = store.submit(make_spec(repo[1], ["ok.py"], labels={"goal": "G1"}))["task_id"]
    assert store.get_task(tid)["state"] == "queued"
    assert worker(store, repos, tmp_path).run_once()
    t = store.get_task(tid)
    assert t["state"] == "succeeded" and t["attempt"] == 1 and t["fence"] == 1
    r = store.get_result(tid)
    assert r["status"] == "succeeded" and r["outcome_class"] == "completed"
    assert r["source"]["resolved_commit"] == repo[1]
    assert r["spec_sha256"] == t["spec_sha256"]
    assert r["runs"][0]["metrics"]["answer"] == 42
    assert r["runs"][0]["user_cpu_seconds"] >= 0 and r["runs"][0]["max_rss_kb"] > 0
    paths = {a["path"] for a in r["artifacts"]}
    assert {"metrics.json", "table.csv", "_taskq_logs/stdout.log"} <= paths
    assert (tmp_path / "results" / tid / "attempt-1-fence-1.json").exists()
    assert store.list_tasks(labels={"goal": "G1"})[0]["task_id"] == tid


def test_benchmark_repetitions_and_summary(store, repos, repo, tmp_path):
    spec = make_spec(repo[1], ["ok.py"], kind="benchmark",
                     benchmark={"warmups": 1, "repetitions": 3})
    spec["command"]["setup"] = [[spec["command"]["argv"][0], "-c", "print('build')"]]
    tid = store.submit(spec)["task_id"]
    worker(store, repos, tmp_path).run_once()
    r = store.get_result(tid)
    assert [x["warmup"] for x in r["runs"]] == [True, False, False, False]
    assert r["summary"]["n"] == 3
    assert r["setup"][0]["exit_code"] == 0 and r["timing"]["setup_wall_seconds"] > 0
    assert {"run-3/metrics.json"} <= {a["path"] for a in r["artifacts"]}


def test_nonzero_exit_is_failed_not_retried(store, repos, repo, tmp_path):
    tid = store.submit(make_spec(repo[1], ["fail.py"]))["task_id"]
    worker(store, repos, tmp_path).run_once()
    r = store.get_result(tid)
    assert r["status"] == "failed" and r["runs"][0]["exit_code"] == 3
    assert store.get_task(tid)["attempt"] == 1


def test_timeout_is_not_completed(store, repos, repo, tmp_path):
    tid = store.submit(make_spec(repo[1], ["sleep.py", "30"],
                                 limits={"timeout_seconds": 0.5}))["task_id"]
    t0 = time.time()
    worker(store, repos, tmp_path).run_once()
    r = store.get_result(tid)
    assert r["status"] == "timeout" and r["outcome_class"] == "not_completed"
    assert time.time() - t0 < 15


def test_unknown_commit_is_infra_error_after_retries(store, repos, repo, tmp_path):
    tid = store.submit(make_spec("b" * 40, ["ok.py"], retry={"max_attempts": 2}))["task_id"]
    w = worker(store, repos, tmp_path)
    w.run_once()
    t = store.get_task(tid)
    assert t["state"] == "queued" and len(t["infra_failures"]) == 1
    w.run_once()
    r = store.get_result(tid)
    assert r["status"] == "infra_error" and r["attempt"] == 2


def test_disallowed_repo_is_infra_error(store, repos, repo, tmp_path):
    spec = make_spec(repo[1], ["ok.py"], retry={"max_attempts": 1})
    spec["source"]["repo"] = "cryptanalysis"
    tid = store.submit(spec)["task_id"]
    worker(store, repos, tmp_path).run_once()
    assert "allowlist" in store.get_result(tid)["error"]


def test_idempotency(store, repo):
    spec = make_spec(repo[1], ["ok.py"], idempotency_key="EXP-1/RUN-1")
    a, b = store.submit(spec), store.submit(spec)
    assert a["task_id"] == b["task_id"] and b["deduplicated"]
    with pytest.raises(Exception, match="different spec"):
        store.submit(make_spec(repo[1], ["fail.py"], idempotency_key="EXP-1/RUN-1"))


def test_cancel_queued(store, repos, repo, tmp_path):
    tid = store.submit(make_spec(repo[1], ["ok.py"]))["task_id"]
    assert store.cancel(tid) == "cancelled"
    assert not worker(store, repos, tmp_path).run_once()
    assert store.get_task(tid)["state"] == "cancelled" and store.get_result(tid) is None


def test_cancel_running(store, repos, repo, tmp_path):
    tid = store.submit(make_spec(repo[1], ["sleep.py", "30"]))["task_id"]
    w = worker(store, repos, tmp_path, lease_seconds=1)
    th = threading.Thread(target=w.run_once)
    th.start()
    while store.get_task(tid)["state"] != "running":
        time.sleep(0.05)
    time.sleep(0.3)
    store.cancel(tid)
    th.join(20)
    assert store.get_result(tid)["status"] == "cancelled"


def test_stale_fence_cannot_write(store, repos, repo, tmp_path):
    """A frozen worker wakes after its lease lapsed and another took over."""
    tid = store.submit(make_spec(repo[1], ["ok.py"]))["task_id"]
    q, msg, _ = store.next_message(["cpu"], "zombie", 60_000, 100)
    fence_a, _ = store.claim(q, msg, tid, "zombie")
    time.sleep(0.3)  # zombie is frozen; lease (0.2s here) lapses
    q2, msg2, tid2 = store.next_message(["cpu"], "w2", 200, 100)
    assert tid2 == tid and msg2 == msg
    fence_b, attempt = store.claim(q2, msg2, tid, "w2")
    assert fence_b == fence_a + 1 and attempt == 2
    assert store.get_task(tid)["infra_failures"][0]["reason"].startswith("worker lost")
    fake = {"schema": "taskq.task-result/v1", "task_id": tid, "attempt": 1,
            "fence": fence_a, "spec_sha256": "0" * 64, "status": "succeeded",
            "outcome_class": "completed", "worker": {"id": "zombie", "hostname": "h"},
            "environment": {}, "source": {"repo": "crypto", "commit": repo[1],
                                          "resolved_commit": repo[1]},
            "timing": {"queued_at": 0, "started_at": 0, "finished_at": 0,
                       "setup_wall_seconds": 0}, "runs": [], "artifacts": []}
    with pytest.raises(FencedError):
        store.complete(q, msg, tid, fence_a, fake)
    assert store.heartbeat(q, msg, tid, "zombie", fence_a) == "fenced"
    store.complete(q2, msg2, tid, fence_b, {**fake, "fence": fence_b, "attempt": 2})
    assert store.get_result(tid)["fence"] == fence_b
    with pytest.raises(FencedError):  # results are write-once
        store.complete(q2, msg2, tid, fence_b, {**fake, "fence": fence_b})


def test_crashed_worker_task_is_reclaimed_and_run(store, repos, repo, tmp_path):
    tid = store.submit(make_spec(repo[1], ["ok.py"]))["task_id"]
    q, msg, _ = store.next_message(["cpu"], "crashed", 60_000, 100)
    store.claim(q, msg, tid, "crashed")  # ... and never heard from again
    time.sleep(0.6)
    w = worker(store, repos, tmp_path, wid="w2", lease_seconds=0.5)
    assert w.run_once()
    r = store.get_result(tid)
    assert r["status"] == "succeeded" and r["attempt"] == 2 and r["worker"]["id"] == "w2"


def test_placement_labels(store, repos, repo, tmp_path):
    tid = store.submit(make_spec(repo[1], ["ok.py"],
                                 placement={"require_labels": {"arch": "arm64"}}))["task_id"]
    assert not worker(store, repos, tmp_path, labels={"arch": "x86_64"}).run_once()
    t = store.get_task(tid)
    assert t["state"] == "queued" and t["attempt"] == 0 and t["declines"] == 1
    assert worker(store, repos, tmp_path, wid="w2", labels={"arch": "arm64"}).run_once()
    assert store.get_result(tid)["status"] == "succeeded"


def test_queue_priority_order(store, repos, repo, tmp_path):
    low = store.submit(make_spec(repo[1], ["ok.py"], queue="cpu"))["task_id"]
    high = store.submit(make_spec(repo[1], ["ok.py"], queue="urgent"))["task_id"]
    w = Worker(store, ["urgent", "cpu"], repos, "w1", block_seconds=0.2)
    w.run_once()
    assert store.get_task(high)["state"] == "succeeded"
    assert store.get_task(low)["state"] == "queued"


def test_cwd_and_stats(store, repos, repo, tmp_path):
    spec = make_spec(repo[1], ["where.py"])
    spec["command"]["cwd"] = "sub"
    tid = store.submit(spec)["task_id"]
    assert store.queue_stats()["cpu"]["stream_length"] == 1
    worker(store, repos, tmp_path).run_once()
    assert store.get_result(tid)["status"] == "succeeded"
    assert store.workers()[0]["id"] == "w1"


def test_shutdown_hands_task_back(store, repos, repo, tmp_path):
    """A pod eviction mid-run requeues the task; it records no result."""
    tid = store.submit(make_spec(repo[1], ["sleep.py", "30"]))["task_id"]
    w = worker(store, repos, tmp_path)
    th = threading.Thread(target=w.run_once)
    th.start()
    while store.get_task(tid)["state"] != "running":
        time.sleep(0.05)
    time.sleep(0.3)
    w.request_shutdown()
    th.join(20)
    t = store.get_task(tid)
    assert t["state"] == "queued" and store.get_result(tid) is None
    assert "shutdown" in t["infra_failures"][0]["reason"]


def test_mcp_tools(store, repos, repo, tmp_path, monkeypatch):
    pytest.importorskip("mcp")
    from taskq import mcp_server
    monkeypatch.setattr(mcp_server, "_store", store)
    sub = mcp_server.submit_command("crypto", repo[1], [__import__("sys").executable, "ok.py"],
                                    repetitions=2, warmups=0, labels={"exp": "EXP-X"})
    assert mcp_server.get_task(sub["task_id"])["spec"]["kind"] == "benchmark"
    worker(store, repos, tmp_path).run_once()
    assert mcp_server.wait_for_task(sub["task_id"], 5)["state"] == "succeeded"
    res = mcp_server.get_result(sub["task_id"], include_runs=False)
    assert "runs" not in res and res["summary"]["n"] == 2
    assert mcp_server.list_tasks(labels={"exp": "EXP-X"})[0]["task_id"] == sub["task_id"]
    assert "task_spec" in mcp_server.get_protocol()


def test_sparse_checkout_and_tree_eviction(store, repo, tmp_path):
    from taskq.execute import RepoCache
    rc = RepoCache(tmp_path / "c", {"crypto": str(repo[0])}, max_trees=1)
    spec = make_spec(repo[1], ["where.py"])
    spec["source"]["sparse_paths"] = ["sub"]
    spec["command"]["cwd"] = "sub"
    tid = store.submit(spec)["task_id"]
    Worker(store, ["cpu"], rc, "w1", block_seconds=0.2).run_once()
    assert store.get_result(tid)["status"] == "succeeded"
    trees = list((tmp_path / "c" / "trees" / "crypto").iterdir())
    assert len(trees) == 1 and (trees[0] / "sub" / "where.py").exists()
    # a second shape at the same commit evicts the first (max_trees=1)
    tid2 = store.submit(make_spec(repo[1], ["ok.py"]))["task_id"]
    Worker(store, ["cpu"], rc, "w1", block_seconds=0.2).run_once()
    assert store.get_result(tid2)["status"] == "succeeded"
    assert [p.name for p in (tmp_path / "c" / "trees" / "crypto").iterdir()] == [repo[1]]

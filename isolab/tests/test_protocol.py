import pytest

from isolab import protocol


def minimal(**over):
    spec = {"schema": "isolab.job/v1", "command": {"argv": ["true"]}}
    spec.update(over)
    return spec


def test_defaults_are_written_out():
    s = protocol.normalize_spec(minimal())
    # normalisation is idempotent: a stored spec re-validates on any worker
    assert protocol.normalize_spec(s) == s
    assert s["pool"] == "default"
    assert s["resources"] == {"cpus": 1, "memory_mb": None, "numa_node": "single", "smt": "isolate", "gpus": 0, "pids": 4096, "scratch_mb": 0}
    assert s["fidelity"]["policy"] == "standard" and s["fidelity"]["min_isolation_tier"] == "C"
    assert s["measure"]["perf_events"] == protocol.DEFAULT_PERF_EVENTS
    assert s["command"]["cwd"] == "." and s["command"]["env"] == {}
    assert s["verify"] is None and s["inputs"] == [] and s["build"] == []


def test_policy_presets_and_overrides():
    s = protocol.normalize_spec(minimal(fidelity={"policy": "strict", "max_other_cpu": 0.5, "require": ["thp_never"]}))
    f = s["fidelity"]
    assert f["min_isolation_tier"] == "A" and f["governor"] == "performance" and f["turbo"] == "off"
    assert f["max_other_cpu"] == 0.5  # explicit wins
    assert "thp_never" in f["require"] and "perf_counters" in f["require"] and "bare_metal" in f["require"]
    b = protocol.normalize_spec(minimal(fidelity={"policy": "best_effort"}))["fidelity"]
    assert b["min_isolation_tier"] == "D" and b["require"] == []


def test_hash_is_order_insensitive_and_default_insensitive():
    a = protocol.normalize_spec(minimal(resources={"cpus": 1}, labels={"b": "2", "a": "1"}))
    b = protocol.normalize_spec(minimal(labels={"a": "1", "b": "2"}))
    assert protocol.spec_sha256(a) == protocol.spec_sha256(b)
    c = protocol.normalize_spec(minimal(resources={"cpus": 2}))
    assert protocol.spec_sha256(c) != protocol.spec_sha256(a)


def test_local_file_inputs_are_refused_at_submission():
    with pytest.raises(protocol.SpecError, match="local_file"):
        protocol.normalize_spec(minimal(inputs=[{"path": "x", "local_file": "/etc/hostname"}]))
    s = protocol.normalize_spec(minimal(inputs=[{"path": "x", "local_file": "/etc/hostname"}]), allow_local_files=True)
    assert s["inputs"][0]["mode"] == "0644"


def test_input_rules():
    with pytest.raises(protocol.SpecError, match="exactly one"):
        protocol.normalize_spec(minimal(inputs=[{"path": "x", "content": "a", "sha256": "0" * 64, "bytes": 1}]))
    with pytest.raises(protocol.SpecError, match="given twice"):
        protocol.normalize_spec(minimal(inputs=[{"path": "x", "content": "a"}, {"path": "./x", "content": "b"}]))
    with pytest.raises(protocol.SpecError, match="relative path"):
        protocol.normalize_spec(minimal(inputs=[{"path": "../x", "content": "a"}]))
    with pytest.raises(protocol.SpecError, match="byte count"):
        protocol.normalize_spec(minimal(inputs=[{"path": "x", "sha256": "0" * 64}]))
    s = protocol.normalize_spec(minimal(inputs=[{"path": "bin/s", "sha256": "0" * 64, "bytes": 3, "mode": "755"}]))
    assert s["inputs"][0]["mode"] == "0755"


def test_cwd_and_backend_rules():
    with pytest.raises(protocol.SpecError, match="relative"):
        protocol.normalize_spec(minimal(command={"argv": ["x"], "cwd": "/abs"}))
    with pytest.raises(protocol.SpecError, match="no image"):
        protocol.normalize_spec(minimal(runtime={"backend": "direct", "image": "x"}))
    with pytest.raises(protocol.SpecError, match="runsc"):
        protocol.normalize_spec(minimal(runtime={"backend": "direct", "oci_runtime": "runsc"}))
    with pytest.raises(protocol.SpecError):
        protocol.normalize_spec(minimal(verify={}))


def test_schema_rejects_unknown_fields():
    with pytest.raises(protocol.SpecError, match="(?i)additional"):
        protocol.normalize_spec(minimal(bogus=1))


def test_job_ids_sort_by_time():
    import time
    a = protocol.new_job_id()
    time.sleep(0.003)
    b = protocol.new_job_id()
    assert a.startswith("J-") and a < b


def test_outcome_class_and_tiers():
    assert protocol.outcome_class("succeeded") == "completed"
    assert protocol.outcome_class("failed") == "completed"
    assert protocol.outcome_class("timeout") == "not_completed"
    assert protocol.tier_at_least("A", "B") and not protocol.tier_at_least("C", "B") and protocol.tier_at_least("D", None)


def test_result_schema_minimal():
    res = {"schema": "isolab.result/v1", "job_id": "J-x", "spec_sha256": "0" * 64, "attempt": 1, "fence": 3,
           "status": "succeeded", "outcome_class": "completed", "error": None, "worker": {"id": "w", "hostname": "h"},
           "host": {}, "placement": {}, "runtime": {}, "fidelity": {"policy": "standard", "grade": "C", "contended": False, "checks": [], "violations": []},
           "timing": {"queued_at": 1, "started_at": 2, "finished_at": 3}, "build": [], "runs": [], "summary": None,
           "artifacts": [], "labels": {}}
    protocol.validate_result(res)
    res["fidelity"]["grade"] = "Z"
    with pytest.raises(Exception):
        protocol.validate_result(res)

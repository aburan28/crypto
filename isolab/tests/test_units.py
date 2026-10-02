import json
import os
import platform
import struct
from pathlib import Path

import pytest

from isolab import perfstat
from isolab.blobs import LocalBlobCache, elf_info, resolve_local_inputs, sha256_bytes
from isolab.fabric import compact_result
from isolab.isolation import Check, grade, max_tier, post_checks, tier_from_mechanisms

PERF_CSV = """# started on Thu
1234.567890,msec,task-clock,1234567890,100.00,0.998,CPUs utilized
42,,context-switches,1234567890,100.00,34.019,/sec
0,,cpu-migrations,1234567890,100.00,0.000,/sec
1500,,page-faults,1234567890,100.00,1.215,K/sec
<not supported>,,cycles,0,100.00,,
<not supported>,,instructions,0,100.00,,
"""

PERF_CSV_HW = """3000000000,,cycles,1000,100.00,2.4,GHz
6000000000,,instructions,1000,100.00,2.00,insn per cycle
100.5,msec,task-clock,1000,100.00,1.0,CPUs utilized
"""


def test_perf_csv_parsing_and_counters():
    parsed = perfstat.parse_csv(PERF_CSV)
    assert parsed["events"]["cycles"]["status"] == "not supported"
    assert parsed["events"]["context-switches"]["value"] == 42.0
    c = perfstat.counters(parsed)
    assert c["task_clock_ms"] == 1234.56789 and c["context_switches"] == 42 and c["cpu_migrations"] == 0
    assert c["instructions"] is None and c["ipc"] is None and c["hardware_counters"] is False
    c2 = perfstat.counters(perfstat.parse_csv(PERF_CSV_HW))
    assert c2["instructions"] == 6000000000 and c2["ipc"] == 2.0 and c2["hardware_counters"] is True
    argv = perfstat.argv_prefix("/usr/bin/perf", ["cycles", "instructions"], [2, 4], "/tmp/o.csv")
    assert argv == ["/usr/bin/perf", "stat", "-x,", "-o", "/tmp/o.csv", "-e", "cycles,instructions", "-a", "-C", "2,4", "--"]


def test_perf_probe_without_perf():
    assert perfstat.probe(None, lambda a, t: (0, "")) == {"available": False, "hw_events": False, "detail": "perf not installed"}
    p = perfstat.probe("/perf", lambda a, t: (0, PERF_CSV))
    assert p["available"] and not p["hw_events"]
    p = perfstat.probe("/perf", lambda a, t: (0, PERF_CSV_HW))
    assert p["available"] and p["hw_events"]


def _elf(machine: int, interp: bytes | None) -> bytes:
    # 64-bit little-endian ELF with one PT_INTERP program header
    phoff = 64
    phentsize, phnum = 56, 1 if interp else 0
    data_off = phoff + phentsize * phnum
    header = b"\x7fELF" + bytes([2, 1, 1, 0]) + b"\0" * 8
    header += struct.pack("<HHIQQQIHHHHHH", 2, machine, 1, 0x400000, phoff, 0, 0, 64, phentsize, phnum, 64, 0, 0)
    ph = b""
    if interp:
        ph = struct.pack("<IIQQQQQQ", 3, 4, data_off, 0, 0, len(interp) + 1, len(interp) + 1, 1)
    return header + ph + (interp + b"\0" if interp else b"")


def test_elf_info(tmp_path):
    p = tmp_path / "dyn"
    p.write_bytes(_elf(0xB7, b"/lib/ld-linux-aarch64.so.1"))
    assert elf_info(p) == {"arch": "aarch64", "bits": 64, "static": False, "interpreter": "/lib/ld-linux-aarch64.so.1"}
    s = tmp_path / "static"
    s.write_bytes(_elf(0x3E, None))
    assert elf_info(s) == {"arch": "x86_64", "bits": 64, "static": True, "interpreter": None}
    t = tmp_path / "text"
    t.write_text("#!/bin/sh\n")
    assert elf_info(t) is None


def test_blob_cache_and_local_inputs(tmp_path):
    cache = LocalBlobCache(tmp_path / "cache")
    d = cache.store_bytes(b"hello")
    assert d == sha256_bytes(b"hello") and cache.has(d)
    cache.path_for(d).write_bytes(b"tampered")
    assert not cache.has(d)
    exe = tmp_path / "tool.sh"
    exe.write_text("#!/bin/sh\necho hi\n")
    exe.chmod(0o755)
    doc = tmp_path / "doc.txt"
    doc.write_text("x")
    spec = {"inputs": [{"path": "bin/tool", "local_file": str(exe)}, {"path": "doc.txt", "local_file": "doc.txt", "mode": "0600"}]}
    uploaded = []
    done = resolve_local_inputs(spec, lambda dg, p: uploaded.append((dg, p.name)), lambda dg: False, base_dir=tmp_path)
    assert [u[1] for u in uploaded] == ["tool.sh", "doc.txt"] and len(done) == 2
    assert spec["inputs"][0]["mode"] == "0755" and spec["inputs"][1]["mode"] == "0600"
    assert "local_file" not in spec["inputs"][0] and spec["inputs"][0]["bytes"] == exe.stat().st_size
    with pytest.raises(FileNotFoundError):
        resolve_local_inputs({"inputs": [{"path": "x", "local_file": "/nonexistent/zzz"}]}, lambda *a: None)


def caps(**over):
    c = {"linux": True, "root": True, "affinity": True, "cpuset_cgroup": True, "cpu_partition": True, "isolated_partition": True,
         "evict": True, "irq_affinity": True, "numa_bind": True, "perf_hw": True}
    c.update(over)
    return c


def test_tiers():
    assert max_tier(caps()) == "A"
    assert max_tier(caps(perf_hw=False)) == "B"
    assert max_tier(caps(perf_hw=False, evict=False)) == "C"
    assert max_tier(caps(linux=False, affinity=False)) == "D"
    m = {"partition_state": "isolated", "evicted": True, "irqs_ok": True, "numa_bound": True, "cgroup": True, "pinned": True}
    assert tier_from_mechanisms(m, caps()) == "A"
    assert tier_from_mechanisms({**m, "partition_state": "root"}, caps()) == "B"
    assert tier_from_mechanisms({**m, "partition_state": "isolated"}, caps(perf_hw=False)) == "B"
    assert tier_from_mechanisms({"pinned": True}, caps()) == "C"
    assert tier_from_mechanisms({}, caps()) == "D"
    assert grade("A", False) == "A" and grade("A", True) == "B" and grade("D", True) == "D" and grade("none", False) == "none"


def summary(**over):
    s = {"n_samples": 5, "wall_s": 5.0, "psi_some_avg10_max": {"cpu": 0.1, "memory": 0.0, "io": 0.3},
         "psi_some_pct": {"cpu": 0.1, "memory": 0.0, "io": 0.3}, "job_cpu_steal_pct": 0.0,
         "job_cpu_busy_pct": 99.0, "irqs_on_job_cpus": 50, "freq_cv_max": 0.001, "throttle_events": 0, "thermal_max_c": 50.0,
         "net_bytes": 0, "contended_samples": 0, "other_cpu_ratio": 0.001, "other_top": {}}
    s.update(over)
    return s


def fid(**over):
    f = {"max_other_cpu": 0.05, "max_psi_some_pct": 1.0, "max_job_cpu_steal_pct": 0.5, "governor": "performance", "require": []}
    f.update(over)
    return f


def test_post_checks_clean_run():
    checks = post_checks(summary(), {"cpu_migrations": 0, "hardware_counters": True}, {"nr_throttled": 0, "oom_kill": 0},
                         {"involuntary_switches": 3, "wall_s": 4.9}, fid(), {"partition_state": "isolated", "reserved_cpus": [2, 3]}, 5.0)
    by = {c.name: c for c in checks}
    assert all(c.status != "fail" for c in checks), [c.to_dict() for c in checks if c.status == "fail"]
    assert by["cpu_migrations"].status == "pass" and by["frequency_stable"].status == "pass"
    assert by["irqs_on_reserved"].status == "pass" and by["irqs_on_reserved"].value == 5.0
    storm = post_checks(summary(irqs_on_job_cpus=50_000), None, None, None, fid(), {"reserved_cpus": [2]}, 5.0)
    assert {c.name: c.status for c in storm}["irqs_on_reserved"] == "fail"
    assert by["launch_overhead_s"].value == pytest.approx(0.1)


def test_post_checks_detect_contention():
    checks = post_checks(summary(other_cpu_ratio=0.3, psi_some_pct={"cpu": 5.0, "memory": 0.0, "io": 0.0},
                                 cg_psi_some_pct={"cpu": 4.0, "memory": 0.0},
                                 job_cpu_steal_pct=3.0, freq_cv_max=0.2, throttle_events=2),
                         {"cpu_migrations": 7, "hardware_counters": False}, {"nr_throttled": 1, "oom_kill": 1}, None,
                         fid(require=["perf_counters"]), {"partition_state": "isolated"}, None)
    failed = {c.name for c in checks if c.status == "fail"}
    assert failed == {"other_cpu", "psi_cpu", "steal", "frequency_stable", "thermal_throttle", "cpu_migrations",
                      "hardware_counters", "cpu_throttled", "oom_kill"}
    # without a partition, migrations are informational
    checks = post_checks(summary(), {"cpu_migrations": 7, "hardware_counters": True}, None, None, fid(), {"partition_state": None}, None)
    assert {c.name: c.status for c in checks}["cpu_migrations"] == "info"
    # host pressure alone does not fail a run whose own cgroup was never stalled
    quiet = post_checks(summary(psi_some_pct={"cpu": 9.0, "memory": 0.0, "io": 0.0}, cg_psi_some_pct={"cpu": 0.0, "memory": 0.0}),
                        None, None, None, fid(), {"reserved_cpus": [2]}, None)
    by = {c.name: c for c in quiet}
    assert by["psi_cpu"].status == "pass" and by["host_psi_cpu"].value == 9.0
    assert post_checks({"n_samples": 0}, None, None, None, fid(), {}, None)[0].status == "unavailable"


def test_compact_result_trims_until_it_fits():
    res = {"host": {"images": ["x" * 5000], "topology": {"cpus": {str(i): {"a": 1} for i in range(500)}}},
           "runs": [{"conditions": {"other_top": {"p": 1}, "freq_mhz": {"0": {}}}, "counters": {"_status": {}}}],
           "fidelity": {"checks": [{"name": "x", "detail": "y" * 3000}]}, "logs": {"stdout_tail": "z" * 9000}}
    out = compact_result(res, 6000)
    assert "images" not in out["host"] and "host.images" in out["truncated"]
    assert len(json.dumps(out)) <= 6000 or "logs" in out["truncated"]


@pytest.mark.skipif(platform.system() != "Linux", reason="reservation mechanics are Linux-only")
def test_reservation_unprivileged_enters_and_exits(tmp_path):
    from isolab.inventory import Paths, cpu_topology
    from isolab.isolation import Reservation, probe_capabilities
    from isolab.planner import plan_cpus
    topo = cpu_topology(Paths())
    online = topo["online"]
    p = plan_cpus(topo, online, 1, "any", "allow")
    c = probe_capabilities()
    fidel = {"policy": "best_effort", "require": [], "max_other_cpu": 1e9, "max_psi_some_pct": 1e9, "max_job_cpu_steal_pct": 100.0,
             "governor": "any", "turbo": "any", "min_isolation_tier": "D"}
    with Reservation("test", p, fidel, c, [x for x in online if x not in p.cpus], None, 100, lock_path=str(tmp_path / "lock")) as r:
        checks = r.pre_checks({"topology": topo, "cpufreq": {"available": False}, "knobs": {}, "virtualization": {}, "cmdline": {}})
        assert any(ch.name == "isolation_tier" for ch in checks)
        settle = r.quiesce(0.3, retries=1, sample_period=0.1)
        assert any(ch.name.startswith("settle") for ch in settle)
        if c["cgroup_writable"]:
            assert r.mechanisms["cgroup"] and r.job_cgroup is not None and r.job_cgroup.exists()
    if c["cgroup_writable"]:
        assert not r.job_cgroup.exists()

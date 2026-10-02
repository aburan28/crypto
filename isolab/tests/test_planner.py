import pytest

from isolab import inventory as inv
from isolab import protocol
from isolab.planner import PlanError, capacity, match_worker, plan_cpus


@pytest.fixture
def topo(synthetic_host):
    return inv.cpu_topology(synthetic_host)


def test_isolate_takes_whole_cores_and_keeps_siblings_idle(topo):
    p = plan_cpus(topo, list(range(2, 16)), 2, "single", "isolate")
    assert p.cpus == [2, 4] and p.idle_siblings == [3, 5] and p.nodes == [0] and p.mems == [0]
    assert p.reserved == [2, 3, 4, 5] and p.cores == [[2, 3], [4, 5]]


def test_allow_fills_cores_first(topo):
    p = plan_cpus(topo, list(range(2, 16)), 3, "single", "allow")
    assert p.cpus == [2, 3, 4] and p.idle_siblings == []


def test_isolate_skips_half_cores(topo):
    # cpu 1's sibling (0) is outside the pool, so core (0,0) is unusable under isolate
    p = plan_cpus(topo, [1, 2, 3], 1, "single", "isolate")
    assert p.cpus == [2] and p.idle_siblings == [3]
    with pytest.raises(PlanError, match="offers at most 1"):
        plan_cpus(topo, [1, 2, 3], 2, "single", "isolate")


def test_numa_node_selection(topo):
    p = plan_cpus(topo, list(range(16)), 2, 1, "isolate")
    assert p.cpus == [8, 10] and p.nodes == [1]
    with pytest.raises(PlanError, match="node 1"):
        plan_cpus(topo, list(range(8)), 1, 1, "isolate")
    assert plan_cpus(topo, list(range(16)), 4, "single", "isolate").cpus == [0, 2, 4, 6]
    with pytest.raises(PlanError, match="any single NUMA node"):
        plan_cpus(topo, list(range(16)), 5, "single", "isolate")
    spread = plan_cpus(topo, list(range(16)), 5, "any", "isolate")
    assert spread.cpus == [0, 2, 4, 6, 8] and spread.nodes == [0, 1] and spread.mems == [0, 1]


def test_smt_off_requires_host_smt_off(topo):
    with pytest.raises(PlanError, match="smt=off requires"):
        plan_cpus(topo, list(range(16)), 1, "single", "off")
    topo2 = {**topo, "smt_control": "off", "summary": {**topo["summary"], "threads_per_core": 1}}
    assert plan_cpus(topo2, list(range(16)), 1, "single", "off").cpus == [0]


def test_busy_cpus_are_skipped(topo):
    p = plan_cpus(topo, list(range(16)), 1, "single", "isolate", busy={0, 1, 2, 3})
    assert p.cpus == [4]


def test_capacity(topo):
    c = capacity(topo, list(range(2, 16)))
    assert c["isolate"] == {"single_node_max": 4, "total": 7, "per_node": {0: 3, 1: 4}}
    assert c["allow"]["single_node_max"] == 8 and c["allow"]["total"] == 14


def worker(topo, **over):
    w = {"id": "w1", "hostname": "h", "pools": ["default"], "labels": {"site": "home"}, "arch": "x86_64",
         "kernel": "6.8.0-1", "cpu": {"flags": ["avx2", "pclmulqdq"]}, "topology": topo, "lab_cpus": list(range(2, 16)),
         "backends": ["podman", "direct"], "default_backend": "podman", "default_image": None,
         "oci_runtimes": ["crun", "runsc"], "images": [{"names": ["localhost/isolab-base:latest"], "digest": "sha256:aa"}],
         "gpus": [], "memory": {"total_kb": 134217728},
         "capabilities": {"perf_hw": True, "cpu_partition": True, "numa_bind": True, "irq_affinity": True, "evict": True, "cpufreq_writable": False},
         "max_tier": "A", "virtualization": {"is_vm": False}, "knobs": {"thp_enabled": "madvise"},
         "cpufreq": {"available": True, "governors": {2: "powersave"}, "turbo": {"control": "/x", "enabled": True}}}
    w.update(over)
    return w


def spec(**over):
    s = {"schema": protocol.SPEC_SCHEMA_ID, "command": {"argv": ["x"]}, "runtime": {"image": "isolab-base"}}
    s.update(over)
    return protocol.normalize_spec(s)


def test_match_worker_accepts_and_explains(topo):
    assert match_worker(spec(), worker(topo)) == []
    r = match_worker(spec(pool="gpu", placement={"labels": {"site": "lab"}, "arch": "aarch64", "cpu_flags": ["avx512f"], "min_kernel": "6.9"}), worker(topo))
    assert any("pool" in x for x in r) and any("label" in x for x in r) and any("arch" in x for x in r)
    assert any("avx512f" in x for x in r) and any("older than 6.9" in x for x in r)
    assert match_worker(spec(placement={"worker": "w2"}), worker(topo)) == ["job is pinned to worker 'w2'"]
    assert match_worker(spec(placement={"worker": "h"}), worker(topo)) == []


def test_match_worker_images_backends_runtimes(topo):
    assert match_worker(spec(runtime={"image": "localhost/isolab-base:latest"}), worker(topo)) == []
    assert match_worker(spec(runtime={"image": "isolab-base:latest"}), worker(topo)) == []
    assert match_worker(spec(runtime={"image": "sha256:aa"}), worker(topo)) == []
    assert "image 'isolab-sage' not present" in match_worker(spec(runtime={"image": "isolab-sage"}), worker(topo))
    assert match_worker(spec(runtime={"image": "isolab-sage"}, placement={"require_images": False}), worker(topo)) == []
    assert any("backend 'docker'" in x for x in match_worker(spec(runtime={"backend": "docker", "image": "isolab-base"}), worker(topo)))
    assert any("runc" in x for x in match_worker(spec(runtime={"oci_runtime": "runc", "image": "isolab-base"}), worker(topo)))
    assert any("no image" in x for x in match_worker(spec(runtime={}), worker(topo)))
    assert match_worker(spec(runtime={"backend": "direct"}), worker(topo)) == []


def test_match_worker_resources_and_fidelity(topo):
    assert any("offers at most 4" in x for x in match_worker(spec(resources={"cpus": 5}), worker(topo)))
    assert any("gpu" in x for x in match_worker(spec(resources={"gpus": 1}), worker(topo)))
    assert any("MiB" in x for x in match_worker(spec(resources={"memory_mb": 70000}), worker(topo)))
    assert match_worker(spec(resources={"memory_mb": 60000}), worker(topo)) == []
    assert any("tier B" in x for x in match_worker(spec(fidelity={"policy": "strict"}), worker(topo, max_tier="B")))
    r = match_worker(spec(fidelity={"policy": "strict"}), worker(topo, virtualization={"is_vm": True}))
    assert any("virtual machine" in x for x in r) and any("governor" in x for x in r) and any("turbo" in x for x in r)
    assert match_worker(spec(fidelity={"policy": "strict", "governor": "any", "turbo": "any"}), worker(topo)) == []
    r = match_worker(spec(fidelity={"require": ["perf_counters", "thp_never", "smt_off"]}), worker(topo, capabilities={"perf_hw": False}))
    assert any("perf_counters" in x for x in r) and any("THP" in x for x in r) and any("SMT" in x for x in r)

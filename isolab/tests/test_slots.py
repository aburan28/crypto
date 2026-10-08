"""Several slots on one host: core splitting, shared host settings, IRQ undo, locks, and concurrent jobs."""
import asyncio
import fcntl
import os
import platform
import sys
import threading
import time

import pytest

from isolab import inventory as inv
from isolab import protocol
from isolab.fabric import Fabric
from isolab.isolation import HostSettings, QuiesceError, Reservation, grade, move_irqs, restore_irqs
from isolab.planner import PlanError, plan_cpus, split_slots
from isolab.worker import SlotRegistry, Worker


def test_split_slots_whole_cores_per_node(synthetic_host):
    topo = inv.cpu_topology(synthetic_host)
    # 16 threads, 8 cores, 2 nodes: four slots of two whole cores, none crossing a node
    slots = split_slots(topo, list(range(16)), 4)
    assert slots == [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15]]
    # uneven: the first slots get the extra core
    slots = split_slots(topo, list(range(2, 16)), 3)
    assert [len(s) for s in slots] == [6, 4, 4]
    assert not set(slots[0]) & set(slots[1]) and sorted(sum(slots, [])) == list(range(2, 16))
    for s in slots:  # never half a core
        assert all((c ^ 1) in s for c in s)
    with pytest.raises(PlanError):
        split_slots(topo, [0, 1], 2)


def test_host_settings_refcount_and_crash(tmp_path):
    knob = tmp_path / "knob"
    knob.write_text("orig")
    hs = HostSettings(str(tmp_path / "state"))
    assert hs.acquire("k", str(knob), "set", f"{os.getpid()}:a") is None
    assert knob.read_text() == "set"
    assert hs.acquire("k", str(knob), "set", f"{os.getpid()}:b") is None
    hs.release("k", f"{os.getpid()}:a")
    assert knob.read_text() == "set"  # b still holds it
    hs.release("k", f"{os.getpid()}:b")
    assert knob.read_text() == "orig"
    # a crashed holder (dead pid) does not pin the value, and the original survives it
    hs.acquire("k", str(knob), "set", "999999999:dead")
    assert knob.read_text() == "set"
    assert hs.acquire("k", str(knob), "set", f"{os.getpid()}:c") is None
    hs.release("k", f"{os.getpid()}:c")
    assert knob.read_text() == "orig"


def _irq_tree(tmp_path, masks):
    proc = tmp_path / "proc"
    for irq, m in masks.items():
        d = proc / "irq" / str(irq)
        d.mkdir(parents=True)
        (d / "smp_affinity_list").write_text(m)
    return proc


def _masks(proc):
    return {d.name: (d / "smp_affinity_list").read_text() for d in sorted((proc / "irq").iterdir())}


@pytest.mark.parametrize("order", ["ab", "ba"])
def test_irq_moves_from_two_slots_undo_in_either_order(tmp_path, order):
    # irq 1 spans both slots and nothing else; irq 2 spans slot a and housekeeping; irq 3 is housekeeping only
    proc = _irq_tree(tmp_path, {1: "2,4", 2: "0,2", 3: "0-1"})
    before = _masks(proc)
    hk = {0, 1}
    a = move_irqs({2, 3}, hk, str(proc))
    b = move_irqs({4, 5}, hk, str(proc))
    assert _masks(proc) == {"1": "0-1", "2": "0", "3": "0-1"}
    for s in (a, b) if order == "ab" else (b, a):
        restore_irqs(s, str(proc))
    assert _masks(proc) == before


def _res(tmp_path, job, exclusive, slot):
    from isolab.planner import Placement
    p = Placement(cpus=[0], idle_siblings=[], nodes=[0], mems=[0], cores=[[0]], smt="allow")
    caps = {"affinity": False, "cpuset_cgroup": False, "evict": False, "irq_affinity": False}
    return Reservation(job, p, {"policy": "standard"}, caps, [], None, 100, lock_path=str(tmp_path / "lock"),
                       lock_wait_s=1.0, slot_lock_path=str(tmp_path / f"lock.{slot}"), exclusive=exclusive)


def test_slot_locks_share_the_host_and_strict_takes_it_whole(tmp_path):
    r1 = _res(tmp_path, "j1", False, "s0")
    r2 = _res(tmp_path, "j2", False, "s1")
    r1._acquire_lock()
    r2._acquire_lock()  # two slots at once
    with pytest.raises(QuiesceError):
        _res(tmp_path, "j3", False, "s0")._acquire_lock()  # the same slot twice
    with pytest.raises(QuiesceError):
        _res(tmp_path, "j4", True, "s1")._acquire_lock()  # the whole host while slots run
    r1.__exit__(None, None, None)
    r2.__exit__(None, None, None)
    strict = _res(tmp_path, "j5", True, "s0")
    strict._acquire_lock()
    assert strict.mechanisms["host_exclusive"] is True
    with pytest.raises(QuiesceError):
        _res(tmp_path, "j6", False, "s1")._acquire_lock()  # no slot starts under a strict job
    strict.__exit__(None, None, None)
    assert not (tmp_path / "lock.drain").exists()


def test_waiting_strict_job_drains_the_slots(tmp_path):
    r1 = _res(tmp_path, "j1", False, "s0")
    r1._acquire_lock()
    strict = _res(tmp_path, "j2", True, "s1")
    strict.lock_wait_s = 10.0
    t = threading.Thread(target=strict._acquire_lock)
    t.start()
    time.sleep(0.3)
    assert (tmp_path / "lock.drain").exists()
    with pytest.raises(QuiesceError, match="whole host is waiting"):
        _res(tmp_path, "j3", False, "s2")._acquire_lock()  # a new slot job backs off
    r1.__exit__(None, None, None)
    t.join(10)
    assert strict._lock is not None and not (tmp_path / "lock.drain").exists()
    strict.__exit__(None, None, None)


def test_grade_caps_shared_host_at_b():
    assert grade("A", False, shared_host=True) == "B"
    assert grade("A", True, shared_host=True) == "B"
    assert grade("C", False, shared_host=True) == "C"
    assert grade("A", False) == "A"


def test_registry_reports_the_other_jobs():
    reg = SlotRegistry()
    reg.enter("a", lambda: {1, 2})
    reg.enter("b", lambda: {3})
    assert reg.others("a") == {"b": {3}}
    reg.leave("b")
    assert reg.others("a") == {}


BUSY = "import time\nt = time.time()\nwhile time.time() - t < 1.5:\n    pass\n"


PRIVILEGED = os.environ.get("ISOLAB_TEST_PRIVILEGED") == "1" and hasattr(os, "geteuid") and os.geteuid() == 0 \
    and platform.system() == "Linux"


@pytest.mark.skipif(platform.system() != "Linux" or (os.cpu_count() or 1) < 3, reason="needs Linux and 3 cpus")
async def test_two_slots_run_concurrently_and_record_each_other(nats_url, fabric, tmp_path):
    await _two_slots(nats_url, fabric, tmp_path, "best_effort")


@pytest.mark.skipif(not PRIVILEGED or (os.cpu_count() or 1) < 3, reason="needs root on Linux with ISOLAB_TEST_PRIVILEGED=1")
async def test_privileged_two_slots_partitions_side_by_side(nats_url, fabric, tmp_path):
    a, b = await _two_slots(nats_url, fabric, tmp_path, "standard")
    for r in (a, b):
        iso = r["placement"]["isolation"]
        assert iso["cgroup"] and iso["evicted"] is True
        assert iso["job_cgroup"].split("/")[-2].startswith("isolab.lab.s")
        # neither job was evicted by the other slot: it had its cpu for the whole loop
        assert r["runs"][0]["cgroup"]["usage_s"] > 1.2
    assert a["placement"]["isolation"]["job_cgroup"] != b["placement"]["isolation"]["job_cgroup"]


async def _two_slots(nats_url, fabric, tmp_path, policy):
    topo = inv.cpu_topology(inv.Paths())
    slots = split_slots(topo, topo["online"][1:], 2)
    host_lab = sorted(sum(slots, []))
    reg = SlotRegistry()
    fabrics, workers, tasks = [], [], []
    for i, cpus in enumerate(slots):
        wf = Fabric(nats_url, namespace=fabric.ns, name=f"test-slot-{i}", lease_s=6.0)
        wf.decline_delay_s = 0.3
        await wf.connect()
        fabrics.append(wf)
        w = Worker(wf, worker_id=f"w.s{i}", pools=["default"], labels={"slot": f"s{i}"}, lab_cpus=cpus,
                   state_dir=tmp_path / "state", backend="direct", lock_path=str(tmp_path / "lock"),
                   block_s=0.5, heartbeat_s=2.0, slot=f"s{i}", host_lab_cpus=host_lab, registry=reg)
        await w.start()
        workers.append(w)
    for w in workers:
        tasks.append(asyncio.create_task(w.run_forever()))
    try:
        ws = await fabric.workers()
        assert sorted(x["id"] for x in ws) == ["w.s0", "w.s1"] and {x["slot"] for x in ws} == {"s0", "s1"}
        ids = []
        for i in range(2):
            s = {"schema": protocol.SPEC_SCHEMA_ID, "runtime": {"backend": "direct"},
                 "inputs": [{"path": "busy.py", "content": BUSY}], "command": {"argv": [sys.executable, "busy.py"]},
                 "placement": {"worker": f"w.s{i}"},
                 "fidelity": {"policy": policy, "settle_s": 0.3}, "limits": {"timeout_s": 60},
                 "measure": {"sample_period_s": 0.2},
                 "resources": {"cpus": 1, "smt": "allow", "numa_node": "any"}}
            ids.append((await fabric.submit(s))["job_id"])
        results = []
        for j in ids:
            rec = await fabric.wait(j, 120)
            assert rec["state"] in protocol.TERMINAL_STATES
            res = await fabric.get_result(j)
            protocol.validate_result(res)
            assert res["status"] == "succeeded", res["error"]
            results.append(res)
        a, b = results
        # the two measured runs overlapped in time
        assert a["timing"]["started_at"] < b["timing"]["finished_at"] and b["timing"]["started_at"] < a["timing"]["finished_at"]
        assert set(a["placement"]["cpus"]).isdisjoint(b["placement"]["cpus"])
        for me, other in ((a, b), (b, a)):
            seen = set(me["runs"][0]["conditions"]["co_tenant_jobs"])
            assert other["job_id"] in seen
            assert me["fidelity"]["shared_host"] is True and me["fidelity"]["grade"] != "A"
            assert me["placement"]["isolation"]["host_exclusive"] is False
            # the other slot's busy loop is a co-tenant, not contention
            co = me["runs"][0]["conditions"]["co_tenant_cpu_s"]
            assert 0.3 < co < me["runs"][0]["conditions"]["wall_s"] + 0.3  # the other job only, never this one
            assert not any(k.startswith("python") and v > 0.3 for k, v in me["runs"][0]["conditions"]["other_top"].items())
        return a, b
    finally:
        for w in workers:
            w.request_shutdown()
        await asyncio.gather(*tasks)
        for wf in fabrics:
            await wf.close()

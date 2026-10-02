"""Which CPUs a job gets, and whether a worker can take a job at all.

Both decisions are pure functions of the spec and the worker's inventory,
so the submitter can predict them from the roster and the worker reproduces
them at claim time.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from . import protocol
from .inventory import cores_of, kernel_version, node_of


class PlanError(ValueError):
    """The requested shape does not fit this host."""


@dataclass
class Placement:
    cpus: list[int]
    idle_siblings: list[int]
    nodes: list[int]
    mems: list[int]
    cores: list[list[int]]
    smt: str
    reason: str = ""
    gpus: list[int] = field(default_factory=list)

    @property
    def reserved(self) -> list[int]:
        return sorted(set(self.cpus) | set(self.idle_siblings))

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["reserved"] = self.reserved
        return d


def _units(topology: dict[str, Any], available: set[int], smt: str) -> list[tuple[int, list[int], list[int]]]:
    """Allocation units as (node, job_cpus, idle_cpus), in deterministic order."""
    units = []
    for key, cpus in sorted(cores_of(topology, available).items()):
        siblings = set(topology["cpus"][cpus[0]]["siblings"])
        whole = siblings <= available
        node = node_of(topology, cpus[0])
        if smt == "isolate":
            if not whole:
                continue  # a sibling the OS may use would share the core
            units.append((node, [cpus[0]], cpus[1:]))
        else:
            for c in cpus:
                units.append((node, [c], []))
    return units


def plan_cpus(topology: dict[str, Any], pool: list[int], want: int,
              numa: str | int = "single", smt: str = "isolate",
              busy: set[int] | frozenset[int] = frozenset()) -> Placement:
    online = set(topology.get("online") or topology["cpus"].keys())
    available = (set(pool) & online) - set(busy)
    if want < 1:
        raise PlanError("cpus must be at least 1")
    if smt == "off":
        ctl = topology.get("smt_control")
        if ctl not in ("off", "forceoff") and topology["summary"].get("threads_per_core", 1) > 1:
            raise PlanError(f"smt=off requires SMT disabled on the host (smt/control is {ctl!r}); "
                            "use smt=isolate to keep siblings idle instead")
    units = _units(topology, available, smt)
    by_node: dict[int, list] = {}
    for u in units:
        by_node.setdefault(u[0], []).append(u)
    nodes = sorted(by_node)
    if numa == "single" or isinstance(numa, int):
        candidates = [n for n in nodes if (numa == "single" or n == numa) and len(by_node[n]) >= want]
        if not candidates:
            best = max((len(v) for v in by_node.values()), default=0)
            where = f"node {numa}" if isinstance(numa, int) else "any single NUMA node"
            raise PlanError(f"wants {want} cpu(s) ({'whole cores' if smt == 'isolate' else 'threads'}) on "
                            f"{where}; this worker offers at most {best} "
                            f"from lab cpus {sorted(available)}")
        chosen = by_node[candidates[0]][:want]
        reason = f"lowest NUMA node with {want} free unit(s)"
    elif numa == "any":
        chosen = []
        for n in nodes:
            chosen.extend(by_node[n][: want - len(chosen)])
            if len(chosen) >= want:
                break
        if len(chosen) < want:
            raise PlanError(f"wants {want} cpu(s); this worker offers {len(chosen)} across all nodes "
                            f"from lab cpus {sorted(available)}")
        reason = "filled across NUMA nodes in id order"
    else:
        raise PlanError(f"numa_node {numa!r} not understood")
    cpus = sorted(c for _, js, _ in chosen for c in js)
    idle = sorted(c for _, _, ids in chosen for c in ids)
    used_nodes = sorted({n for n, _, _ in chosen})
    cores = [sorted(js + ids) for _, js, ids in chosen]
    return Placement(cpus=cpus, idle_siblings=idle, nodes=used_nodes, mems=used_nodes,
                     cores=cores, smt=smt, reason=reason)


def capacity(topology: dict[str, Any], pool: list[int]) -> dict[str, Any]:
    """How big a job this pool can take, per SMT mode, on one node and overall."""
    out: dict[str, Any] = {}
    for smt in ("isolate", "allow"):
        units = _units(topology, set(pool) & set(topology.get("online") or pool), smt)
        per_node: dict[int, int] = {}
        for n, _, _ in units:
            per_node[n] = per_node.get(n, 0) + 1
        out[smt] = {"single_node_max": max(per_node.values(), default=0), "total": len(units),
                    "per_node": per_node}
    return out


def _image_matches(wanted: str, have: list[dict[str, Any]]) -> bool:
    if not wanted:
        return False
    w = wanted
    if w.startswith("sha256:"):
        return any((img.get("digest") or "") == w or (img.get("id") or "").startswith(w[7:]) for img in have)
    if ":" not in w.rsplit("/", 1)[-1]:
        w = w + ":latest"
    for img in have:
        for name in img.get("names") or []:
            if name == w or name.endswith("/" + w):
                return True
            if w.endswith("/" + name):
                return True
    return False


def match_worker(spec: dict[str, Any], inv: dict[str, Any]) -> list[str]:
    """Reasons this worker cannot take the normalised spec; empty means it can."""
    reasons: list[str] = []
    place, res, fid, rt = spec["placement"], spec["resources"], spec["fidelity"], spec["runtime"]
    if spec["pool"] not in inv.get("pools", []):
        reasons.append(f"not in pool {spec['pool']!r}")
    if place["worker"] and place["worker"] not in (inv.get("id"), inv.get("hostname")):
        reasons.append(f"job is pinned to worker {place['worker']!r}")
    for k, v in (place.get("labels") or {}).items():
        if inv.get("labels", {}).get(k) != v:
            reasons.append(f"lacks label {k}={v}")
    if place["arch"] and place["arch"] != inv.get("arch"):
        reasons.append(f"arch {inv.get('arch')} is not {place['arch']}")
    flags = set(inv.get("cpu", {}).get("flags") or [])
    missing = [f for f in place["cpu_flags"] if f not in flags]
    if missing:
        reasons.append(f"cpu lacks flags {missing}")
    if place["min_kernel"]:
        want = tuple(int(x) for x in place["min_kernel"].split("."))
        if kernel_version(inv.get("kernel")) < want:
            reasons.append(f"kernel {inv.get('kernel')} is older than {place['min_kernel']}")
    backends = inv.get("backends") or []
    backend = rt["backend"]
    if backend == "auto":
        backend = inv.get("default_backend") or (backends[0] if backends else None)
        if backend is None:
            reasons.append("no execution backend available")
    elif backend not in backends:
        reasons.append(f"backend {backend!r} not available (has {backends})")
    if rt["oci_runtime"] != "auto" and rt["oci_runtime"] not in (inv.get("oci_runtimes") or []):
        reasons.append(f"OCI runtime {rt['oci_runtime']!r} not available (has {inv.get('oci_runtimes')})")
    if backend and backend != "direct":
        image = rt["image"] or inv.get("default_image")
        if not image:
            reasons.append("no image given and the worker has no default image")
        elif place["require_images"] and not _image_matches(image, inv.get("images") or []):
            reasons.append(f"image {image!r} not present")
    if res["gpus"] > len(inv.get("gpus") or []):
        reasons.append(f"wants {res['gpus']} gpu(s); worker has {len(inv.get('gpus') or [])}")
    topo = inv.get("topology")
    if topo:
        try:
            # inventories travel as JSON, so cpu keys may be strings
            topo = _intkeys(topo)
            p = plan_cpus(topo, inv.get("lab_cpus") or [], res["cpus"], res["numa_node"], res["smt"])
            if res["memory_mb"]:
                if res["numa_node"] != "any" and len(p.nodes) == 1:
                    node_kb = topo["nodes"][p.nodes[0]].get("mem_total_kb")
                else:
                    node_kb = inv.get("memory", {}).get("total_kb")
                if node_kb and res["memory_mb"] * 1024 > node_kb:
                    reasons.append(f"wants {res['memory_mb']} MiB; node(s) {p.nodes} have {node_kb // 1024} MiB")
        except PlanError as err:
            reasons.append(str(err))
    caps = inv.get("capabilities") or {}
    have_tier = inv.get("max_tier", "D")
    if not protocol.tier_at_least(have_tier, fid.get("min_isolation_tier")):
        reasons.append(f"can reach isolation tier {have_tier}, job needs {fid['min_isolation_tier']}")
    need_caps = {"perf_counters": "perf_hw", "cpu_partition": "cpu_partition", "numa_bind": "numa_bind",
                 "irq_moved": "irq_affinity", "evicted": "evict", "frequency_pinned": "cpufreq_writable"}
    for req in fid.get("require") or []:
        if req in need_caps and not caps.get(need_caps[req]):
            reasons.append(f"requires {req}, which this worker cannot provide")
        elif req == "bare_metal" and inv.get("virtualization", {}).get("is_vm"):
            reasons.append("requires bare metal; this worker is a virtual machine")
        elif req == "thp_never" and inv.get("knobs", {}).get("thp_enabled") != "never":
            reasons.append(f"requires THP never; host has {inv.get('knobs', {}).get('thp_enabled')}")
        elif req == "smt_off" and topo and topo.get("smt_control") not in ("off", "forceoff") \
                and topo.get("summary", {}).get("threads_per_core", 1) > 1:
            reasons.append("requires SMT off; host has SMT enabled")
    cf = inv.get("cpufreq") or {}
    if fid.get("governor") == "performance" and cf.get("available") and not caps.get("cpufreq_writable"):
        govs = set((cf.get("governors") or {}).values())
        if govs and govs != {"performance"}:
            reasons.append(f"governor must be performance; host runs {sorted(govs)} and the worker cannot change it")
    if fid.get("turbo") == "off" and cf.get("turbo", {}).get("control") and not caps.get("cpufreq_writable"):
        if cf["turbo"].get("enabled"):
            reasons.append("turbo must be off; it is on and the worker cannot change it")
    return reasons


def _intkeys(topo: dict[str, Any]) -> dict[str, Any]:
    out = dict(topo)
    out["cpus"] = {int(k): v for k, v in topo["cpus"].items()}
    out["nodes"] = {int(k): v for k, v in topo["nodes"].items()}
    return out

from isolab import inventory as inv
from isolab import samplers


def test_cpulist_round_trip():
    assert inv.parse_cpulist("0-3,8,10-11") == [0, 1, 2, 3, 8, 10, 11]
    assert inv.format_cpulist([0, 1, 2, 3, 8, 10, 11]) == "0-3,8,10-11"
    assert inv.parse_cpulist("") == [] and inv.format_cpulist([]) == ""
    assert inv.format_cpulist({5}) == "5"


def test_topology(synthetic_host):
    t = inv.cpu_topology(synthetic_host)
    assert t["summary"] == {"sockets": 2, "cores": 8, "threads_per_core": 2, "numa_nodes": 2}
    assert t["cpus"][9] == {"core_id": 0, "package_id": 1, "node": 1, "siblings": [8, 9], "online": True}
    assert t["nodes"][1]["cpus"] == list(range(8, 16)) and t["nodes"][1]["mem_total_kb"] == 67108864
    assert t["smt_control"] == "on" and t["caches"][2]["shared_cpu_list"] == "0-7"
    assert inv.cores_of(t, [0, 1, 2, 9]) == {(0, 0): [0, 1], (0, 1): [2], (1, 0): [9]}


def test_cpuinfo_cmdline_knobs(synthetic_host):
    ci = inv.cpuinfo(synthetic_host)
    assert ci["model"].startswith("Intel") and "avx512f" in ci["flags"] and ci["microcode"] == "0x2b000603"
    cmd = inv.kernel_cmdline(synthetic_host)
    assert cmd["isolcpus"] == [2, 3, 4, 5, 6, 7] and cmd["isolcpus_flags"] == ["managed_irq", "domain"]
    assert cmd["nohz_full"] == [2, 3, 4, 5, 6, 7] and cmd["mitigations"] == "auto"
    k = inv.kernel_knobs(synthetic_host)
    assert k["thp_enabled"] == "madvise" and k["aslr"] == 2 and k["numa_balancing"] == 1 and k["nmi_watchdog"] == 1
    cf = inv.cpufreq(synthetic_host, [0, 1])
    assert cf["available"] and cf["governors"] == {0: "performance", 1: "powersave"}
    assert cf["turbo"]["enabled"] is True and cf["turbo"]["control"].endswith("no_turbo")
    assert inv.memory(synthetic_host)["total_kb"] == 134217728
    assert inv.psi_available(synthetic_host)
    assert inv.cgroups(synthetic_host)["controllers"] == ["cpuset", "cpu", "io", "memory", "pids"]


def test_samplers_readers(synthetic_host):
    psi = samplers.read_psi("io", synthetic_host)
    assert psi["some"]["avg10"] == 0.5 and psi["full"]["total"] == 900
    st = samplers.read_cpu_stat(synthetic_host)
    assert st[3]["steal"] == 2 and st[0]["idle"] == 100
    irqs = samplers.read_interrupts(synthetic_host)
    assert irqs[0] == 101 and irqs[15] == 116
    assert samplers.read_softirqs(synthetic_host)[4] == 7
    assert samplers.read_freq([1], synthetic_host) == {1: 2401000}
    assert samplers.read_netdev(synthetic_host) == {"rx_bytes": 50000, "tx_bytes": 20000}
    assert samplers.read_diskstats(synthetic_host) == {"read_sectors": 8000, "write_sectors": 4000}
    assert samplers.read_loadavg(synthetic_host) == [0.1, 0.2, 0.3]


def test_images_parsing():
    rows = '[{"Names":["localhost/isolab-base:latest"],"Id":"abc","Digest":"sha256:ff","Size":1}]'
    out = inv.images(lambda argv, t: (0, rows))
    assert out[0]["names"] == ["localhost/isolab-base:latest"] and out[0]["digest"] == "sha256:ff"
    lines = '{"Repository":"python","Tag":"3.12","ID":"deadbeef","Size":"1MB"}\n'
    out = inv.images(lambda argv, t: (0, lines), backend="docker")
    assert out[0]["names"] == ["python:3.12"]


def test_gpus_parsing():
    text = "0, NVIDIA RTX PRO 6000, GPU-abc, 98304, 580.65, Enabled, Default, 2617, 14001, 00000000:01:00.0\n"
    out = inv.gpus(lambda argv, t: (0, text), which=lambda n: "/usr/bin/nvidia-smi")
    assert out[0]["name"] == "NVIDIA RTX PRO 6000" and out[0]["memory_mb"] == 98304 and out[0]["persistence_mode"] == "Enabled"
    assert inv.gpus(lambda a, t: (0, ""), which=lambda n: None) == []

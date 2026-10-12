import asyncio
import os
import shutil
import socket
import subprocess
import time
from pathlib import Path

import pytest

from isolab.fabric import Fabric


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="session")
def nats_url(tmp_path_factory):
    """A JetStream server: ISOLAB_TEST_NATS_URL if set, else one started from ISOLAB_TEST_NATS_SERVER / PATH."""
    if os.environ.get("ISOLAB_TEST_NATS_URL"):
        yield os.environ["ISOLAB_TEST_NATS_URL"]
        return
    exe = os.environ.get("ISOLAB_TEST_NATS_SERVER") or shutil.which("nats-server")
    if not exe:
        pytest.skip("no nats-server: set ISOLAB_TEST_NATS_URL or ISOLAB_TEST_NATS_SERVER")
    port = _free_port()
    store = tmp_path_factory.mktemp("js")
    proc = subprocess.Popen([exe, "-js", "-a", "127.0.0.1", "-p", str(port), "-sd", str(store)],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    url = f"nats://127.0.0.1:{port}"
    for _ in range(100):
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.2):
                break
        except OSError:
            time.sleep(0.05)
    yield url
    proc.terminate()
    proc.wait()


_counter = {"n": 0}


@pytest.fixture
async def fabric(nats_url, request):
    _counter["n"] += 1
    ns = f"T{abs(hash(request.node.nodeid)) % 10**7}{_counter['n']}"
    f = Fabric(nats_url, namespace=ns, name="test-client", lease_s=5.0)
    f.decline_delay_s = 0.3
    await f.connect()
    yield f
    await f.close()


@pytest.fixture
async def fabric2(nats_url, fabric):
    """A second connection to the same namespace (another worker or client)."""
    f = Fabric(nats_url, namespace=fabric.ns, name="test-client-2", lease_s=5.0)
    f.decline_delay_s = 0.3
    await f.connect()
    yield f
    await f.close()


def write(p: Path, text: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


@pytest.fixture
def synthetic_host(tmp_path):
    """A /sys and /proc tree for a 2-socket, 8-core, 16-thread, 2-NUMA-node machine.

    cpu c: package c//8, node c//8, core (c//2)%4, siblings (c&~1, c|1).
    """
    sys_ = tmp_path / "sys"
    proc = tmp_path / "proc"
    cpu = sys_ / "devices/system/cpu"
    write(cpu / "present", "0-15\n")
    write(cpu / "online", "0-15\n")
    write(cpu / "smt/control", "on\n")
    for c in range(16):
        t = cpu / f"cpu{c}/topology"
        write(t / "core_id", f"{(c // 2) % 4}\n")
        write(t / "physical_package_id", f"{c // 8}\n")
        write(t / "thread_siblings_list", f"{c & ~1}-{c | 1}\n")
        write(cpu / f"cpu{c}/online", "1\n")
        cf = cpu / f"cpu{c}/cpufreq"
        write(cf / "scaling_governor", "powersave\n" if c else "performance\n")
        write(cf / "scaling_driver", "intel_pstate\n")
        write(cf / "scaling_cur_freq", f"{2400000 + c * 1000}\n")
        write(cf / "cpuinfo_min_freq", "800000\n")
        write(cf / "cpuinfo_max_freq", "3800000\n")
        write(cf / "scaling_available_governors", "performance powersave\n")
    write(cpu / "intel_pstate/no_turbo", "0\n")
    for lvl, typ, size, shared in ((1, "Data", "48K", "0-1"), (2, "Unified", "2048K", "0-1"), (3, "Unified", "32768K", "0-7")):
        d = cpu / f"cpu0/cache/index{lvl}"
        write(d / "level", f"{lvl}\n"); write(d / "type", f"{typ}\n"); write(d / "size", f"{size}\n"); write(d / "shared_cpu_list", f"{shared}\n")
    for n in range(2):
        nd = sys_ / f"devices/system/node/node{n}"
        write(nd / "cpulist", f"{n * 8}-{n * 8 + 7}\n")
        write(nd / "meminfo", f"Node {n} MemTotal:       67108864 kB\nNode {n} MemFree:        50000000 kB\n")
        write(nd / "distance", "10 21\n" if n == 0 else "21 10\n")
    write(sys_ / "kernel/mm/transparent_hugepage/enabled", "always [madvise] never\n")
    write(sys_ / "kernel/mm/transparent_hugepage/defrag", "always defer defer+madvise [madvise] never\n")
    write(sys_ / "kernel/mm/ksm/run", "0\n")
    write(sys_ / "fs/cgroup/cgroup.controllers", "cpuset cpu io memory pids\n")
    write(sys_ / "fs/cgroup/cgroup.subtree_control", "cpuset cpu memory pids\n")
    write(sys_ / "class/dmi/id/sys_vendor", "Supermicro\n")
    write(sys_ / "class/dmi/id/product_name", "X13SAE\n")
    write(proc / "cpuinfo", "processor\t: 0\nvendor_id\t: GenuineIntel\nmodel name\t: Intel(R) Xeon(R) w5-2455X\nmicrocode\t: 0x2b000603\n"
                            "flags\t\t: fpu sse2 avx2 pclmulqdq avx512f\n\nprocessor\t: 1\nmodel name\t: Intel(R) Xeon(R) w5-2455X\n\n")
    write(proc / "cmdline", "BOOT_IMAGE=/vmlinuz root=/dev/sda1 isolcpus=managed_irq,domain,2-7 nohz_full=2-7 rcu_nocbs=2-7 mitigations=auto\n")
    write(proc / "meminfo", "MemTotal:       134217728 kB\nMemAvailable:   100000000 kB\nSwapTotal:             0 kB\nSwapFree:              0 kB\n")
    write(proc / "sys/kernel/randomize_va_space", "2\n")
    write(proc / "sys/kernel/perf_event_paranoid", "-1\n")
    write(proc / "sys/kernel/numa_balancing", "1\n")
    write(proc / "sys/kernel/nmi_watchdog", "1\n")
    write(proc / "sys/kernel/timer_migration", "1\n")
    write(proc / "sys/vm/swappiness", "60\n")
    write(proc / "pressure/cpu", "some avg10=0.00 avg60=0.10 avg300=0.05 total=123456\nfull avg10=0.00 avg60=0.00 avg300=0.00 total=0\n")
    write(proc / "pressure/memory", "some avg10=0.00 avg60=0.00 avg300=0.00 total=10\nfull avg10=0.00 avg60=0.00 avg300=0.00 total=5\n")
    write(proc / "pressure/io", "some avg10=0.50 avg60=0.10 avg300=0.05 total=999\nfull avg10=0.40 avg60=0.00 avg300=0.00 total=900\n")
    stat = "cpu  100 0 50 1000 5 2 3 0 0 0\n" + "".join(f"cpu{c} 10 0 5 100 0 1 1 {2 if c == 3 else 0} 0 0\n" for c in range(16)) + "intr 5000\nctxt 1000\n"
    write(proc / "stat", stat)
    header = "           " + "".join(f"CPU{c:<10}" for c in range(16)) + "\n"
    write(proc / "interrupts", header + " 24:" + "".join(f"{c + 1:>11}" for c in range(16)) + "  IR-PCI-MSI  eth0\n"
          + "LOC:" + "".join(f"{100:>11}" for c in range(16)) + "  Local timer interrupts\n"
          + "ERR:          0\n")
    write(proc / "softirqs", header + "   TIMER:" + "".join(f"{7:>11}" for c in range(16)) + "\n")
    write(proc / "loadavg", "0.10 0.20 0.30 1/200 12345\n")
    write(proc / "net/dev", "Inter-|   Receive                                                |  Transmit\n face |bytes    packets errs drop fifo frame compressed multicast|bytes    packets errs drop fifo colls carrier compressed\n"
                            "    lo:    1000      10    0    0    0     0          0         0     1000      10    0    0    0     0       0          0\n"
                            "  eth0:   50000     100    0    0    0     0          0         0    20000      50    0    0    0     0       0          0\n")
    write(proc / "diskstats", " 259       0 nvme0n1 100 0 8000 10 50 0 4000 5 0 0 0 0 0 0 0 0 0\n 259       1 nvme0n1p1 1 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n")
    from isolab.inventory import Paths
    return Paths(sys=str(sys_), proc=str(proc), etc=str(tmp_path / "etc"))

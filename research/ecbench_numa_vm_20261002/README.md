# ecbench on a two-node Linux guest: pinning and NUMA binding

**Result: on a two-node Linux kernel, `ecbench` pins every measured
child to its reserved CPU and binds its memory to that CPU's node.**
The policy reads back as `bind:<node>`, and every anonymous page the
solve touched sits on that node. Reaching that took three attempts. The
second found a defect the single-node CI runner cannot see (the NUMA
read-back checked the cpuset, not the policy), and it was fixed in
`494aaa8fa`. All three attempts are kept here.

This is a functional test of the isolation machinery, not a
performance measurement. The guest is emulated (QEMU TCG), so its wall
times mean nothing and its remote-memory latency is not real. What it
can show, and does, is that the kernel accepted every placement
`ecbench` asked for and that `ecbench` reads each one back correctly.

## Setup

| | |
|---|---|
| guest | `qemu-system-x86_64 -M q35 -accel tcg,thread=multi -cpu max,vendor=GenuineIntel -smp 8,sockets=2,cores=2,threads=2`, 2 × 512 MiB NUMA nodes (CPUs 0–3, 4–7); full command in [`vm/run-vm.sh`](vm/run-vm.sh) |
| kernel | Ubuntu `6.8.0-142-generic` (`noble-server-cloudimg-amd64-vmlinuz-generic`, SHA-256 `cd5fcfd260b91782637b7b4e221a48e656358f6eef549602540e62380b6f2f2c`), cmdline `psi=1 isolcpus=6,7` |
| userland | an initramfs holding busybox 1.35.0 x86_64-musl (SHA-256 `6e123e7f3202a8c1e9b1f94d8941580a25135382b99e8d3e34fb858bba311348`), `ecbench`, [`vm/spec.json`](vm/spec.json) and [`vm/init`](vm/init). Init runs everything as root, then powers off |
| binary | `ecbench` from `494aaa8fa`, cross-built for `x86_64-unknown-linux-musl` (static-pie) with LLVM clang and `ld.lld`, SHA-256 `3ff9988dd054d0222f17afea61cae00b944b466c0358b5e582a86eff922921d1` |
| spec | a prime and a Koblitz curve (`icv1-fp20-t727-cd198a38`, `icv1-f2m17-tm101-00378d4e`), 2 targets each, arms `rho.negation`, `bsgs.negation`, `kangaroo.vow` and an A/A `rho.negation`, 3 rounds after 1 warm-up: 64 executions per session |

The init script runs three sessions:

1. `--cpus auto`, which should take the highest whole core away from
   CPU 0: CPUs 6 and 7, on node 1, both in `isolcpus`.
2. `--cpus 2,3`, a whole core on node 0.
3. `--cpus 2`, a lone SMT sibling, which must be refused.

## Results (the fixed binary)

[`sessions/auto`](sessions/auto) and [`sessions/node0`](sessions/node0),
console in [`vm/console-summary.txt`](vm/console-summary.txt):

| session | run CPU | node | policy read back | anonymous pages | affinity read back | CPU at start, end | levels | verified |
|---|---:|---:|---|---|---|---|---|---:|
| `--cpus auto` | 6 | 1 | `bind:1` | 100 % on node 1 | `[6]` | 6, 6 | L2 ×51, L1 ×13 | 64 / 64 |
| `--cpus 2,3` | 2 | 0 | `bind:0` | 100 % on node 0 | `[2]` | 2, 2 | L2 ×52, L1 ×12 | 64 / 64 |
| `--cpus 2` | — | — | — | — | — | — | refused: "CPU 2's SMT sibling 3 is not reserved" (exit 2) | — |

- The guest showed the topology the command asked for: `node0 cpus=0-3`,
  `node1 cpus=4-7`, siblings `0-1`, `2-3`, `6-7`, isolated `6-7`, and PSI
  for cpu, io and memory.
- Every run that stopped at L1 did so on the run-queue delay check alone
  (0.8–10.8 % of the solve against a 0.5 % limit), which is emulator
  contention.
  Every run failed L3 on the governor, turbo and bare-metal checks, as a
  VM must. The `isolcpus` check passed on CPU 6 and failed on CPU 2,
  which is correct.
- **Audits.** Each session audits OK in the guest and again on macOS
  arm64, with 12 of 12 replays identical (`audit-*-on-macos.json`). The
  file hashes agree between the two audits, so the serial transfer was
  byte-exact, and counts recorded on x86-64 Linux reproduce exactly on
  arm64 macOS.
- **Caveat.** TCG reports zero steal by construction, so these L2 grades
  certify the guest kernel's view and nothing about the host underneath.

## The attempts that failed, kept

1. **Alpine `virt` kernel 6.6**
   ([`vm/console-summary-alpine-no-numa.txt`](vm/console-summary-alpine-no-numa.txt)).
   That kernel is built without `CONFIG_NUMA` (`/sys/devices/system/node`
   absent), and `-cpu max` reported no SMT. It was not a NUMA test, and
   the lone-sibling run was rightly accepted because CPU 2 had no
   sibling. It also exposed a bug: the host capsule said "6 logical" on
   an 8-CPU guest, because `available_parallelism` follows the process's
   affinity and `isolcpus=6,7` narrowed it. The class counted CPUs by how
   ecbench was launched. Fixed in `494aaa8fa`: online CPUs are counted
   from the topology, and the process's affinity is recorded outside the
   class. No session from this attempt is kept, only the console.
2. **Ubuntu kernel, before the fix**
   ([`sessions/before-fix-auto`](sessions/before-fix-auto),
   [`sessions/before-fix-node0`](sessions/before-fix-node0),
   [`vm/console-summary-before-fix.txt`](vm/console-summary-before-fix.txt)).
   Pinning was exact, but every run was graded **L0** with "memory allowed
   on Some("0-1"), not bound to node 1". The read-back compared
   `Mems_allowed_list` with the node. That list is the cpuset's
   permission, not the task's policy, so it reads `0-1` whatever policy
   the task has. Fixed in `494aaa8fa`: the child reads the policy with
   `get_mempolicy` and counts anonymous pages per node in
   `/proc/self/numa_maps`. These sessions were made by a binary built from
   a working tree between `3da9dc08a` and `494aaa8fa` (SHA-256
   `a947d480943a0b85850abe5218c6f0c1462aeb711ed278aea4872e8f8c6ddf47`,
   recorded in each record). Their counts still replay exactly.

A third defect showed in the same consoles: the guest's class id changed
between boots (`ECBENV1he85c57801186`, `ECBENV1hb8c5c33506f0`) because
version 1 of the class hashed each node's `MemTotal`, which moved by
44 MB between boots. Version 2 (`ECBENV2h…`) hashes the topology's shape
and memory in whole GiB. Version-1 ids, these included, still recompute
under their own definition.

## Reproduce

Cross-build the binary (static, no C toolchain for Linux beyond LLVM):

```bash
CC_x86_64_unknown_linux_musl=/opt/homebrew/opt/llvm/bin/clang CFLAGS_x86_64_unknown_linux_musl=--target=x86_64-unknown-linux-musl AR_x86_64_unknown_linux_musl=/opt/homebrew/opt/llvm/bin/llvm-ar RUSTFLAGS="-C linker=/opt/homebrew/opt/lld@21/bin/ld.lld -C linker-flavor=ld.lld -C link-self-contained=yes" cargo build --release --target x86_64-unknown-linux-musl --bin ecbench
```

Then boot the guest with the kernel and busybox above:

```bash
research/ecbench_numa_vm_20261002/vm/run-vm.sh target/x86_64-unknown-linux-musl/release/ecbench vmlinuz-generic busybox
```

Audit the committed sessions anywhere:

```bash
./target/release/ecbench verify --dir research/ecbench_numa_vm_20261002/sessions/auto --replay 12 --exit-code
```

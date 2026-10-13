# Baseline build and host manifest

Command outputs captured before candidate implementation.

```text
rustc 1.99.0 (b940084d7 2026-09-28)
binary: rustc
commit-hash: b940084d7eb6a299eb4bfeb8e34901bc051e7ac4
commit-date: 2026-09-28
host: x86_64-unknown-linux-gnu
release: 1.99.0
LLVM version: 23.1.1
99cf6d805d50d04a0ca578f4346d8399cb05e8ea2608e4b91ed87f0b28dff826  /tmp/crypto-baseline-target/release/ecbench
72ac315b7dd26a2bfa97707a17a3a9103991b65260accc4ca73c2b4a53838b3b  /tmp/crypto-baseline-target/release/examples/koblitz_sixsum_extract_probe
759978250c1ce53d7b87bd8a1c2304781fcc4fc4
Linux 38726d0dfae6 6.18.44 #1 SMP Sat Sep 26 20:02:31 UTC 2026 x86_64 GNU/Linux
Architecture:                            x86_64
CPU op-mode(s):                          32-bit, 64-bit
Address sizes:                           46 bits physical, 57 bits virtual
Byte Order:                              Little Endian
CPU(s):                                  5
On-line CPU(s) list:                     0-4
Vendor ID:                               GenuineIntel
Model name:                              INTEL(R) XEON(R) PLATINUM 8573C
CPU family:                              6
Model:                                   207
Thread(s) per core:                      1
Core(s) per socket:                      5
Socket(s):                               1
Stepping:                                2
BogoMIPS:                                4600.00
Flags:                                   fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ss ht syscall nx pdpe1gb rdtscp lm constant_tsc rep_good nopl xtopology nonstop_tsc cpuid tsc_known_freq pni pclmulqdq ssse3 fma cx16 pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand hypervisor lahf_lm abm 3dnowprefetch cpuid_fault fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm avx512f avx512dq rdseed adx smap avx512ifma clflushopt clwb avx512cd sha_ni avx512bw avx512vl xsaveopt xsavec xgetbv1 xsaves avx_vnni avx512_bf16 arat avx512vbmi umip waitpkg avx512_vbmi2 gfni vaes vpclmulqdq avx512_vnni avx512_bitalg avx512_vpopcntdq rdpid cldemote movdiri movdir64b fsrm serialize tsxldtrk avx512_fp16 arch_capabilities
Hypervisor vendor:                       KVM
Virtualization type:                     full
L1d cache:                               144 KiB (3 instances)
L1i cache:                               96 KiB (3 instances)
L2 cache:                                6 MiB (3 instances)
L3 cache:                                260 MiB (1 instance)
NUMA node(s):                            1
NUMA node0 CPU(s):                       0-4
Vulnerability Gather data sampling:      Not affected
Vulnerability Ghostwrite:                Not affected
Vulnerability Indirect target selection: Vulnerable
Vulnerability Itlb multihit:             Not affected
Vulnerability L1tf:                      Not affected
Vulnerability Mds:                       Not affected
               total        used        free      shared  buff/cache   available
Mem:            17Gi       1.9Gi        15Gi       1.2Gi       1.7Gi        15Gi
Swap:             0B          0B          0B
```

Initial build attempt with the inherited PATH failed before compilation:

```text
/bin/bash: line 1: cargo: command not found
```

A temporary stable toolchain was installed outside the repository under `/tmp`; the binaries above were built from the unmodified baseline revision with `CARGO_TARGET_DIR=/tmp/crypto-baseline-target`.


#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)

= N83 retained factor-base run report

*Run ID:* `run-20261016T2352Z`  \
*Date:* 2026-10-16  \
*Protocol/source:* `research/n83_prebuilt_factor_base_20261016/PROTOCOL.md`,
commit `2d3a7f3723925311e74a29ff1cdc7c74a0a97f50`.  \
*Host:* shared macOS arm64, L0.  \
*Workload:* 54 native compressed objects per panel, 83,328,324 bytes;
both exact degree-83 curves, three policies, three seeds and K=64/256/600.

#image("panel-times.svg", width: 100%)

All three native panel manifests report 54/54 complete, with internal
construction intervals of 59.517962583 s, 74.022112708 s and
169.554772250 s. The corresponding outer real intervals were 60.07 s,
74.73 s and 172.78 s. Every object hash and point-set hash matched the
historical S3 panel exactly.

The first three replay attempts failed because temporary objects were
cleared across a chat continuation. A separate workspace reconstruction
recreated byte-identical objects. Three independent processes then replayed
the original manifests with the generic multi-limb arithmetic path:
46.66 s, 52.15 s and 31.43 s; all returned PASS for 54 objects.

The `/usr/bin/time -l` wrapper itself exited 1 after each completed panel
because `sysctl kern.clockrate` was unavailable in the sandbox. Peak RSS is
missing. Replay and construction do different work, and shared-host wall
spread is large. Relation yield, rank, target work, complete cost, rho and
an equal-output external comparator were not measured.

*Raw and replay receipts:* `research/n83_prebuilt_factor_base_20261016/runs/`.
The main `RESULT.md` gives the exact checksum binding and interpretation.

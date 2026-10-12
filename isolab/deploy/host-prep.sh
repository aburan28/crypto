#!/usr/bin/env bash
# Apply the runtime-tunable host settings for high-fidelity runs and print the boot-time ones.
# Usage: sudo deploy/host-prep.sh --cpus 4-15
set -euo pipefail
CPUS=""; while [ $# -gt 0 ]; do case "$1" in --cpus) CPUS="$2"; shift;; esac; shift; done
[ "$(id -u)" = 0 ] || { echo "run as root" >&2; exit 1; }
isolab doctor ${CPUS:+--cpus $CPUS} --apply
cat >/etc/sysctl.d/90-isolab.conf <<SYS
kernel.nmi_watchdog = 0
kernel.numa_balancing = 0
kernel.timer_migration = 0
kernel.perf_event_paranoid = -1
vm.stat_interval = 120
SYS
sysctl -q -p /etc/sysctl.d/90-isolab.conf || true
echo
echo "Boot-time settings (edit /etc/default/grub GRUB_CMDLINE_LINUX, then update-grub and reboot):"
isolab doctor ${CPUS:+--cpus $CPUS} --json | python3 -c 'import json,sys; [print("  "+i["fix"]) for i in json.load(sys.stdin) if i.get("reboot") and i["status"]!="ok"]'

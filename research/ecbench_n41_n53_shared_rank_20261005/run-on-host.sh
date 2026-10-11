#!/usr/bin/env bash
# L2 host run for the frozen n37 and n41/n53 shared-rank ecbench specs.
#
# Preconditions (see L2_RUNBOOK.md): an already provisioned Linux bare-metal
# host, Ubuntu 24.04, run as root, with outbound network for rustup, crates.io
# and GitHub. Launching the host is the repository owner's decision and is not
# performed by this script.
#
# Usage: COMMIT=<sha> ./run-on-host.sh [WORKDIR]
set -euo pipefail

: "${COMMIT:?set COMMIT to the exact crypto commit to measure}"
REPO_URL="${REPO_URL:-https://github.com/aburan28/crypto.git}"
WORK="${1:-/root/ecbench-l2}"
N37_SPEC="research/ecbench_n37_native_online_wall_20261004/SPEC.json"
N41_N53_SPEC="research/ecbench_n41_n53_shared_rank_20261005/SPEC.json"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

if [ "$(id -u)" -ne 0 ]; then
  echo "run as root: core reservation and thread eviction need it (ecbench README section 6)" >&2
  exit 2
fi
if ! grep -qi 'ubuntu' /etc/os-release || ! grep -q 'VERSION_ID="24.04"' /etc/os-release; then
  echo "warning: this runbook was written for Ubuntu 24.04; continuing on $(. /etc/os-release; echo "$PRETTY_NAME")" >&2
fi

mkdir -p "$WORK"
cd "$WORK"

# 1. Toolchain and build dependencies.
export DEBIAN_FRONTEND=noninteractive
apt-get update -y
apt-get install -y --no-install-recommends build-essential git curl ca-certificates pkg-config \
  linux-tools-common "linux-tools-$(uname -r)" cpufrequtils sqlite3 || true
if ! command -v cargo >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
fi
# shellcheck disable=SC1091
. "$HOME/.cargo/env"

# 2. Exact source.
if [ ! -d crypto ]; then
  git clone --no-checkout "$REPO_URL" crypto
fi
git -C crypto fetch origin "$COMMIT"
git -C crypto checkout --detach "$COMMIT"
test "$(git -C crypto rev-parse HEAD)" = "$COMMIT"

# 3. Build the harness once; the binary hash is the identity every record carries.
( cd crypto && cargo build --release --bin ecbench )
BIN="$WORK/crypto/target/release/ecbench"
sha256sum "$BIN" | tee "$WORK/ecbench.sha256"
rustc --version | tee "$WORK/rustc.version"

# 4. Quiet the host. Every step is best effort and recorded; ecbench's own
#    preflight and per-run grading decide the level each run earns.
{
  systemctl stop unattended-upgrades.service apt-daily.timer apt-daily-upgrade.timer \
    apt-daily.service apt-daily-upgrade.service fwupd-refresh.timer motd-news.timer \
    man-db.timer e2scrub_all.timer 2>&1 || true
  systemctl disable --now unattended-upgrades.service 2>&1 || true
  systemctl stop snapd.service snapd.socket snapd.seeded.service 2>&1 || true
  systemctl mask snapd.service snapd.socket 2>&1 || true
  for g in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
    [ -w "$g" ] && echo performance > "$g" || true
  done
  if [ -w /sys/devices/system/cpu/intel_pstate/no_turbo ]; then echo 1 > /sys/devices/system/cpu/intel_pstate/no_turbo || true; fi
  if [ -w /sys/devices/system/cpu/cpufreq/boost ]; then echo 0 > /sys/devices/system/cpu/cpufreq/boost || true; fi
  swapoff -a || true
  sysctl -w kernel.numa_balancing=0 || true
} 2>&1 | tee "$WORK/quiet-$STAMP.log"

# 5. Host facts beside the sessions.
"$BIN" host | tee "$WORK/HOST-$STAMP.txt"
{ uname -a; lscpu; cat /proc/cmdline; cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null; } > "$WORK/PROVENANCE-$STAMP.txt" 2>&1 || true

# 6. Run the frozen specs. No --allow-busy: a busy preflight refuses to start,
#    which is the point. --cpus auto reserves a whole core away from CPU 0.
SESSIONS="$WORK/sessions-$STAMP"
mkdir -p "$SESSIONS"
cd "$WORK/crypto"
"$BIN" plan --spec "$N37_SPEC" --json > "$SESSIONS/n37-plan.json"
"$BIN" plan --spec "$N41_N53_SPEC" --json > "$SESSIONS/n41-n53-plan.json"
"$BIN" run --spec "$N37_SPEC" --out "$SESSIONS/n37_linux_l2_01" --cpus auto --wait 2>&1 | tee "$SESSIONS/n37-run.log"
"$BIN" run --spec "$N41_N53_SPEC" --out "$SESSIONS/n41_n53_linux_l2_01" --cpus auto --wait 2>&1 | tee "$SESSIONS/n41-n53-run.log"

# 7. Audit with every deterministic measured run replayed.
"$BIN" verify --dir "$SESSIONS/n37_linux_l2_01" --replay-all --exit-code --out "$SESSIONS/n37-AUDIT.json" | tail -n 5
"$BIN" verify --dir "$SESSIONS/n41_n53_linux_l2_01" --replay-all --exit-code --out "$SESSIONS/n41-n53-AUDIT.json" | tail -n 5
"$BIN" table --dir "$SESSIONS/n37_linux_l2_01" > "$SESSIONS/n37-TABLE.txt"
"$BIN" table --dir "$SESSIONS/n41_n53_linux_l2_01" > "$SESSIONS/n41-n53-TABLE.txt"

# 8. Seal.
cp "$WORK/ecbench.sha256" "$WORK/rustc.version" "$WORK/HOST-$STAMP.txt" "$WORK/PROVENANCE-$STAMP.txt" "$WORK/quiet-$STAMP.log" "$SESSIONS/"
( cd "$SESSIONS" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS )
tar -C "$WORK" -czf "$WORK/ecbench-l2-sessions-$STAMP.tar.gz" "sessions-$STAMP"
sha256sum "$WORK/ecbench-l2-sessions-$STAMP.tar.gz"
echo "sessions sealed in $WORK/ecbench-l2-sessions-$STAMP.tar.gz; copy it off the host before terminating"

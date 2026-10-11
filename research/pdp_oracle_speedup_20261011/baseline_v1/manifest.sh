#!/usr/bin/env bash
# Host + toolchain + build manifest for PS-0.  Usage: manifest.sh BINARY... > manifest.json
set -euo pipefail
repo="$(cd "$(dirname "$0")/../../.." && pwd)"
esc() { jq -Rs 'rtrimstr("\n")'; }
flags=$(grep -m1 '^flags' /proc/cpuinfo | tr ' ' '\n' | grep -E '^(avx|avx2|avx512[a-z_0-9]*|bmi1|bmi2|pclmulqdq|vpclmulqdq|gfni|vaes|adx|fma|aes|sha_ni)$' | sort | tr '\n' ' ')
echo '{'
echo "  \"date_utc\": \"$(date -u +%FT%TZ)\","
echo "  \"commit\": \"$(git -C "$repo" rev-parse HEAD)\","
echo "  \"worktree_dirty\": $(git -C "$repo" status --porcelain -- src Cargo.toml Cargo.lock examples | grep -q . && echo true || echo false),"
echo "  \"cpu_model\": $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | esc),"
echo "  \"cpu_flags_of_interest\": $(echo "$flags" | esc),"
echo "  \"vcpus\": $(nproc),"
echo "  \"lscpu_caches\": $(lscpu | grep -E 'L1d|L1i|L2|L3|Hypervisor|Virtualization type' | esc),"
echo "  \"cpufreq_control\": $([ -d /sys/devices/system/cpu/cpu0/cpufreq ] && echo true || echo false),"
echo "  \"cgroup_cpu_max\": $( (cat /sys/fs/cgroup/cpu.max 2>/dev/null || echo none) | esc),"
echo "  \"memory\": $(grep -m1 MemTotal /proc/meminfo | esc),"
echo "  \"kernel\": $(uname -srv | esc),"
echo "  \"rustc\": $(rustc -Vv | esc),"
echo "  \"cargo\": $(cargo --version | esc),"
echo "  \"valgrind\": $(valgrind --version | esc),"
echo "  \"cargo_profile_release\": $(awk '/^\[profile.release\]/{f=1;next} /^\[/{f=0} f' "$repo/Cargo.toml" | esc),"
echo "  \"gpu\": $( (nvidia-smi -L 2>/dev/null || echo none) | esc),"
echo '  "binaries": {'
first=1
for b in "$@"; do
  [ $first = 1 ] || echo ','
  first=0
  printf '    %s: %s' "$(printf '%s' "$(basename "$b")" | esc)" "$(sha256sum "$b" | cut -d' ' -f1 | esc)"
done
echo
echo '  }'
echo '}'

#!/usr/bin/env bash
# PROTOCOL.md part B, on the AWS host: build ecbench from the shipped tree,
# record the host, run eight one-target m = 83 sessions side by side (one
# physical core each), then stop the instance.  Thin orchestration only;
# every number comes from ecbench's own sealed records.
#
#   run-on-host.sh probe            # build, provenance, 10^6-step throughput probe
#   run-on-host.sh start HOURS      # schedule a stop after HOURS, start the 8 sessions
set -euo pipefail
cd "$(dirname "$0")/../.."
R=research/ecbench_m83_rho_20261005
OUT="$R/sessions"
mkdir -p "$OUT" "$R/host"

case "${1:-}" in
probe)
  cargo build --release --bin ecbench
  {
    echo "commit $(cat COMMIT 2>/dev/null || echo unknown)"
    rustc -Vv
    uname -a
    lscpu
    free -g
    sha256sum target/release/ecbench
    sha256sum "$R"/specs/*.json
  } > "$R/host/PROVENANCE.txt"
  ./target/release/ecbench host > "$R/host/HOST.json"
  # Throughput probe: one target, the same method, capped at 10^6 steps
  # (step_cap_factor 0); it exhausts by design and is not a result.
  sed -e 's/"dp_bits": "12"/"dp_bits": "12", "step_cap_factor": "0"/' \
      "$R/specs/target-1.json" > /tmp/probe.json
  rm -rf /tmp/probe-out
  taskset -c 1 ./target/release/ecbench run --spec /tmp/probe.json --out /tmp/probe-out \
      --cpus inherit --lock /tmp/ecbench-probe.lock
  jq -r '"probe wall_s=\(.time.solve_wall_ns/1e9)"' /tmp/probe-out/records.jsonl \
      | tee "$R/host/PROBE.txt"
  ;;
start)
  hours="${2:?hours before the instance stops itself}"
  sudo shutdown -h "+$((hours * 60))" "ecbench m83: hard stop after ${hours} h"
  for i in 1 2 3 4 5 6 7 8; do
    core=$((i - 1))
    nohup taskset -c "$core" ./target/release/ecbench run \
        --spec "$R/specs/target-$i.json" --out "$OUT/m83-target-$i" \
        --cpus inherit --lock "/tmp/ecbench-$i.lock" \
        > "$R/host/run-target-$i.log" 2>&1 &
    echo $! > "$R/host/pid-$i"
  done
  # Stop the instance (shutdown behaviour: stop) once every session ends.
  nohup bash -c "
    for i in 1 2 3 4 5 6 7 8; do
      while kill -0 \$(cat $R/host/pid-\$i) 2>/dev/null; do sleep 60; done
    done
    date -u > $R/host/DONE
    sync
    sudo shutdown -h now 'ecbench m83: all sessions ended'
  " > "$R/host/watcher.log" 2>&1 &
  echo "started; watcher pid $!"
  ;;
*)
  echo "usage: $0 probe | start HOURS" >&2
  exit 2
  ;;
esac

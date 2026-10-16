#!/bin/bash
# Rendered by allocate.py; the boot-time watchdog enforces the approved cap.
set -Eeuo pipefail

RUN_ROOT=/var/lib/n83-tier1
EVIDENCE="$RUN_ROOT/evidence"
BUCKET=__BUCKET__
SOURCE_COMMIT=__COMMIT__
RUN_TOKEN=__RUN_TOKEN__
S3_KEY="n83-tier1-online-20261016/$RUN_TOKEN.tar.gz"
mkdir -p "$EVIDENCE/bin"
exec > >(tee -a /var/log/n83-tier1-bootstrap.log "$EVIDENCE/bootstrap.log") 2>&1

# This timer counts from EC2 pending time, before package installation.
if ! REMAINING_SECONDS=$(python3 - <<'PY'
from datetime import datetime, timezone
import json
import urllib.request

request = urllib.request.Request(
    "http://169.254.169.254/latest/api/token", data=b"", method="PUT",
    headers={"X-aws-ec2-metadata-token-ttl-seconds": "21600"})
with urllib.request.urlopen(request, timeout=3) as response:
    token = response.read().decode()
request = urllib.request.Request(
    "http://169.254.169.254/latest/dynamic/instance-identity/document",
    headers={"X-aws-ec2-metadata-token": token})
with urllib.request.urlopen(request, timeout=3) as response:
    identity = json.load(response)
pending = datetime.fromisoformat(identity["pendingTime"].replace("Z", "+00:00"))
elapsed = int((datetime.now(timezone.utc) - pending).total_seconds())
print(max(0, 165 * 60 - elapsed))
PY
); then
  /usr/sbin/shutdown -h now
  exit 1
fi
if [ "$REMAINING_SECONDS" -lt 60 ]; then
  /usr/sbin/shutdown -h now
  exit 1
fi
if ! systemd-run --unit=n83-tier1-hard-stop \
  --on-active="${REMAINING_SECONDS}s" \
  /usr/sbin/shutdown -h now; then
  /usr/sbin/shutdown -h now
  exit 1
fi
systemctl list-timers n83-tier1-hard-stop.timer > "$EVIDENCE/watchdog.txt"
printf '%s\n' "$REMAINING_SECONDS" > "$EVIDENCE/watchdog-remaining-seconds.txt"
date -u +%FT%TZ > "$EVIDENCE/instance-start-utc.txt"

archive_and_stop() {
  local status=$?
  trap - EXIT
  set +e
  date -u +%FT%TZ > "$EVIDENCE/bootstrap-finished-utc.txt"
  printf '%s\n' "$status" > "$EVIDENCE/bootstrap-exit-status.txt"
  cp /var/log/cloud-init-output.log "$EVIDENCE/cloud-init-output.log" 2>/dev/null
  tar -C "$RUN_ROOT" -czf "$RUN_ROOT/$RUN_TOKEN.tar.gz" evidence
  local archive_status=$?
  if [ "$archive_status" -eq 0 ]; then
    for attempt in 1 2 3; do
      if aws s3api put-object --bucket "$BUCKET" --key "$S3_KEY" \
        --body "$RUN_ROOT/$RUN_TOKEN.tar.gz" \
        --server-side-encryption AES256; then
        echo "ARCHIVED s3://$BUCKET/$S3_KEY"
        /usr/sbin/shutdown -h now
        return 0
      fi
      echo "artifact upload attempt $attempt failed"
      sleep 10
    done
  fi
  echo "ARCHIVE_UPLOAD_FAILED: instance remains available until watchdog termination"
  return 1
}
trap archive_and_stop EXIT

cat /etc/os-release > "$EVIDENCE/os-release.txt"
uname -a > "$EVIDENCE/uname.txt"
cat /proc/cpuinfo > "$EVIDENCE/cpuinfo.txt"
cat /proc/meminfo > "$EVIDENCE/meminfo.txt"

dnf install -y git gcc gcc-c++ make pkgconf-pkg-config openssl-devel curl
if ! command -v aws >/dev/null 2>&1; then
  dnf install -y awscli-2 || dnf install -y awscli
fi
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
  sh -s -- -y --profile minimal --default-toolchain 1.99.0
export PATH="/root/.cargo/bin:$PATH"
rustc --version > "$EVIDENCE/rustc-version.txt"
cargo --version > "$EVIDENCE/cargo-version.txt"
aws --version > "$EVIDENCE/aws-version.txt" 2>&1

git clone --depth 1 --branch codex/n83-tier1-online-20261016 \
  https://github.com/aburan28/crypto.git "$RUN_ROOT/src"
cd "$RUN_ROOT/src"
test "$(git rev-parse HEAD)" = "$SOURCE_COMMIT"
git rev-parse HEAD > "$EVIDENCE/source-commit.txt"
python3 research/n83_tier1_online_20261016/run.py --fixture-preflight > \
  "$EVIDENCE/fixture-preflight.json"

export CARGO_TARGET_DIR="$RUN_ROOT/cargo-target"
export CARGO_BUILD_JOBS=8
export RUST_TEST_THREADS=1
cargo test --release --lib --locked > "$EVIDENCE/libtest.log" 2>&1
cargo test --release --locked --example koblitz_orbit_dlp_fast > \
  "$EVIDENCE/ic-example-tests.log" 2>&1
cargo test --release --locked --example n83_rho_optimization_eval > \
  "$EVIDENCE/rho-example-tests.log" 2>&1
cargo build --release --locked --example koblitz_orbit_dlp_fast > \
  "$EVIDENCE/ic-build.log" 2>&1
cargo build --release --locked --example n83_rho_optimization_eval > \
  "$EVIDENCE/rho-build.log" 2>&1
cp "$CARGO_TARGET_DIR/release/examples/koblitz_orbit_dlp_fast" \
  "$EVIDENCE/bin/koblitz_orbit_dlp_fast"
cp "$CARGO_TARGET_DIR/release/examples/n83_rho_optimization_eval" \
  "$EVIDENCE/bin/n83_rho_optimization_eval"

export PYTHONDONTWRITEBYTECODE=1
timeout 6600 python3 research/n83_tier1_online_20261016/run.py \
  --binary "$EVIDENCE/bin/koblitz_orbit_dlp_fast" \
  --rho-binary "$EVIDENCE/bin/n83_rho_optimization_eval" \
  --output "$EVIDENCE/results" > "$EVIDENCE/run.log" 2>&1
timeout 1200 python3 research/n83_tier1_online_20261016/verify.py \
  "$EVIDENCE/results" > "$EVIDENCE/verify.log" 2>&1

#!/bin/bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
exec 9>field_launcher.lock
flock -n 9 || { echo 'Field replay launcher is already active'; exit 3; }
if pgrep -f "^$root/field_worker $root " >/dev/null; then echo 'Field workers already active'; exit 3; fi
if [[ -e FIELD_STOP ]]; then echo 'FIELD_STOP exists'; exit 3; fi
sha256sum -c FIELD_FROZEN.sha256 --quiet
mkdir -p logs
logroot=$(mktemp -d "$root/logs/field_launch_$(date -u +%Y%m%dT%H%M%SZ)_XXXXXX")
exec node "$root/source/field_controller.mjs" "$root" "$logroot"

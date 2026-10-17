#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
ulimit -v 12582912
exec taskset -c 0 /workspace/native-deps/cm-followups-v1 campaign

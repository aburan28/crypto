#!/usr/bin/env bash
# Install an isolab hub and/or worker on a Linux host. Run as root.
#
#   sudo deploy/install.sh --hub --worker --cpus 4-15 [--pool P]... [--label K=V]... [--gvisor]
#   sudo deploy/install.sh --worker --cpus 2-7 --hub-url nats://TOKEN@first-host:4222
#   sudo deploy/install.sh --hub --cluster-routes nats://first-host:6222
set -euo pipefail
PREFIX=/opt/isolab
HUB=0 WORKER=0 CPUS="" HUB_URL="" GVISOR=1 POOLS=() LABELS=() ROUTES="" DEFAULT_IMAGE="localhost/isolab-base:latest" SRC=""
while [ $# -gt 0 ]; do
  case "$1" in
    --hub) HUB=1;;
    --worker) WORKER=1;;
    --cpus) CPUS="$2"; shift;;
    --pool) POOLS+=("$2"); shift;;
    --label) LABELS+=("$2"); shift;;
    --hub-url) HUB_URL="$2"; shift;;
    --cluster-routes) ROUTES="$2"; shift;;
    --default-image) DEFAULT_IMAGE="$2"; shift;;
    --no-gvisor) GVISOR=0;;
    --gvisor) GVISOR=1;;
    --prefix) PREFIX="$2"; shift;;
    --source) SRC="$2"; shift;;
    -h|--help) sed -n 2,8p "$0"; exit 0;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac; shift
done
[ "$(id -u)" = 0 ] || { echo "run as root" >&2; exit 1; }
[ $HUB = 1 ] || [ $WORKER = 1 ] || { echo "give --hub and/or --worker" >&2; exit 2; }
HERE="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${SRC:-$HERE}"

echo "== packages"
if command -v apt-get >/dev/null; then
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -q
  apt-get install -y -q python3-venv python3-pip gcc libc6-dev git curl ca-certificates numactl hwloc-nox \
    podman crun uidmap slirp4netns fuse-overlayfs >/dev/null
  apt-get install -y -q "linux-tools-$(uname -r)" >/dev/null 2>&1 || apt-get install -y -q linux-tools-generic >/dev/null 2>&1 || echo "perf: install linux-tools for $(uname -r) by hand"
elif command -v dnf >/dev/null; then
  dnf install -y -q python3 gcc glibc-static git curl numactl hwloc podman crun perf >/dev/null
fi

echo "== python environment in $PREFIX"
mkdir -p "$PREFIX" /etc/isolab /var/lib/isolab/hub
python3 -m venv "$PREFIX"
"$PREFIX/bin/pip" install -q --upgrade pip
"$PREFIX/bin/pip" install -q "$SRC[mcp]"
ln -sf "$PREFIX/bin/isolab" /usr/local/bin/isolab

echo "== nats-server"
"$PREFIX/bin/isolab" nats-download --dest "$PREFIX/bin" >/dev/null
"$PREFIX/bin/nats-server" --version

if [ $GVISOR = 1 ] && ! command -v runsc >/dev/null; then
  echo "== gVisor"
  ARCH=$(uname -m); REL=$(curl -fsSL https://api.github.com/repos/google/gvisor/releases/latest | sed -n 's/.*"tag_name": "\([^"]*\)".*/\1/p')
  # the release tarball holds runsc, the containerd shim and a gvisor-bin/ directory of sidecars
  # (gvisor_sentry and friends) that runsc looks for next to itself; install all of it
  TMP=$(mktemp -d); curl -fsSL -o "$TMP/gv.tbz" "https://github.com/google/gvisor/releases/download/$REL/gvisor-$ARCH.tar.bz2" \
    && tar -xjf "$TMP/gv.tbz" -C "$TMP" && rm -f "$TMP/gv.tbz" && cp -a "$TMP"/. /usr/local/bin/ && chown -R root:root /usr/local/bin/runsc /usr/local/bin/containerd-shim-runsc-v1 /usr/local/bin/gvisor-bin && rm -rf "$TMP" \
    && runsc --version | head -1 || echo "gVisor install failed; continue without it"
fi

echo "== launcher and calibration kernel"
"$PREFIX/bin/isolab" launcher-build --dest /var/lib/isolab/bin

id isolab >/dev/null 2>&1 || useradd --system --home /var/lib/isolab --shell /usr/sbin/nologin isolab
chown -R isolab:isolab /var/lib/isolab/hub

TOKEN=""
if [ $HUB = 1 ]; then
  echo "== hub"
  if [ -f /etc/isolab/hub.env ]; then TOKEN=$(sed -n 's/^ISOLAB_TOKEN=//p' /etc/isolab/hub.env); fi
  TOKEN="${TOKEN:-$(python3 -c 'import secrets; print(secrets.token_urlsafe(32))')}"
  EXTRA=""; [ -n "$ROUTES" ] && EXTRA="--cluster-name isolab --routes $ROUTES"
  cat > /etc/isolab/hub.env <<ENV
ISOLAB_HUB_LISTEN=0.0.0.0:4222
ISOLAB_HUB_STORE=/var/lib/isolab/hub
ISOLAB_HUB_NAME=$(hostname -s)
ISOLAB_TOKEN=$TOKEN
ISOLAB_HUB_EXTRA=$EXTRA
ENV
  chmod 600 /etc/isolab/hub.env
  install -m 0644 "$HERE/deploy/isolab-hub.service" /etc/systemd/system/
  systemctl daemon-reload; systemctl enable --now isolab-hub.service
  HUB_URL="${HUB_URL:-nats://$TOKEN@127.0.0.1:4222}"
fi

if [ $WORKER = 1 ]; then
  echo "== worker"
  [ -n "$CPUS" ] || { echo "--cpus is required for a worker (e.g. --cpus 4-15)" >&2; exit 2; }
  [ -n "$HUB_URL" ] || { echo "--hub-url is required when this host runs no hub" >&2; exit 2; }
  EXTRA=""; for p in "${POOLS[@]:-}"; do [ -n "$p" ] && EXTRA="$EXTRA --pool $p"; done
  for l in "${LABELS[@]:-}"; do [ -n "$l" ] && EXTRA="$EXTRA --label $l"; done
  EXTRA="$EXTRA --default-image $DEFAULT_IMAGE"
  cat > /etc/isolab/worker.env <<ENV
ISOLAB_URL=$HUB_URL
ISOLAB_CPUS=$CPUS
ISOLAB_WORKER_EXTRA=$EXTRA
ENV
  chmod 600 /etc/isolab/worker.env
  install -m 0644 "$HERE/deploy/isolab-worker.service" /etc/systemd/system/
  systemctl daemon-reload; systemctl enable --now isolab-worker.service
fi

echo
echo "done. clients connect with:"
echo "  export ISOLAB_URL=${HUB_URL/127.0.0.1/$(hostname -I 2>/dev/null | awk '{print $1}')}"
echo "next:  sudo isolab doctor --cpus ${CPUS:-<lab cpus>}    then    sudo isolab images build base"

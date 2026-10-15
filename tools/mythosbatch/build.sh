#!/bin/sh
# Build static mythos-batch binaries into dist/ (no cgo, no runtime dependencies).
#   ./build.sh            all platforms
#   ./build.sh native     just this machine -> ./mythos-batch
set -eu
cd "$(dirname "$0")"
flags="-trimpath -ldflags=-s -ldflags=-w"
if [ "${1:-}" = native ]; then
  CGO_ENABLED=0 go build $flags -o mythos-batch .
  echo "built ./mythos-batch"
  exit 0
fi
mkdir -p dist
for target in linux/amd64 linux/arm64 darwin/arm64 darwin/amd64 windows/amd64; do
  os=${target%/*}; arch=${target#*/}
  ext=""; [ "$os" = windows ] && ext=".exe"
  out="dist/mythos-batch-$os-$arch$ext"
  CGO_ENABLED=0 GOOS=$os GOARCH=$arch go build $flags -o "$out" .
  echo "built $out"
done

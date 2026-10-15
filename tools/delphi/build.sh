#!/bin/sh
# Build static delphi binaries (no cgo, no runtime dependencies).
#   ./build.sh            all platforms into dist/
#   ./build.sh native     this machine -> ./delphi
set -eu
cd "$(dirname "$0")"
if [ "${1:-}" = native ]; then
  CGO_ENABLED=0 go build -trimpath -ldflags="-s -w" -o delphi ./cmd/delphi
  echo "built ./delphi"
  exit 0
fi
mkdir -p dist
for target in linux/amd64 linux/arm64 darwin/arm64 darwin/amd64 windows/amd64; do
  os=${target%/*}; arch=${target#*/}
  ext=""; [ "$os" = windows ] && ext=".exe"
  out="dist/delphi-$os-$arch$ext"
  CGO_ENABLED=0 GOOS=$os GOARCH=$arch go build -trimpath -ldflags="-s -w" -o "$out" ./cmd/delphi
  echo "built $out"
done
(cd dist && sha256sum delphi-* > SHA256SUMS 2>/dev/null || shasum -a 256 delphi-* > SHA256SUMS)

#!/bin/sh
# Offline Linux arm64 build; compilation is preparation, not a pilot trial.
set -eu
checkout=$(git rev-parse --show-toplevel)
target=/private/tmp/n83-structured-capacity-linux-target
registry=${CARGO_HOME:-$HOME/.cargo}/registry
image=sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31
mkdir -p "$target"
docker run --rm --pull never --network none --cpus 4 \
  --memory 8192m --memory-swap 8192m --pids-limit 256 \
  --cap-drop ALL --security-opt no-new-privileges \
  --user "$(id -u):$(id -g)" --tmpfs /cargo:rw,mode=1777 \
  --mount "type=bind,src=$registry,dst=/cargo/registry,readonly" \
  --mount "type=bind,src=$checkout,dst=/src,readonly" \
  --mount "type=bind,src=$target,dst=/target" \
  --env CARGO_HOME=/cargo --env CARGO_TARGET_DIR=/target \
  --workdir /src "$image" \
  cargo build --release --offline --locked --no-default-features \
  --example koblitz_n83_structured_orbit_capacity

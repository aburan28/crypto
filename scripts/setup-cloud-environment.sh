#!/usr/bin/env bash
# Paste this into the cloud environment setup field as:
#   bash /workspace/scripts/setup-cloud-environment.sh
# Installs system math/build tools plus the repository's Conductor and cairn
# binaries. Do not put tokens or credentials in this file.
set -Eeuo pipefail

repo_dir="${CRYPTO_REPO_DIR:-/workspace}"
bin_dir="$HOME/.local/bin"
mkdir -p "$bin_dir"
export PATH="$bin_dir:$HOME/.cargo/bin:$PATH"

run_root() {
  if [[ "$(id -u)" -eq 0 ]]; then
    "$@"
  else
    if ! command -v sudo >/dev/null 2>&1; then
      echo "sudo is required to install system packages" >&2
      return 1
    fi
    sudo "$@"
  fi
}

if ! command -v apt-get >/dev/null 2>&1; then
  echo "This setup script requires a Debian or Ubuntu image with apt-get" >&2
  exit 1
fi

packages=(
  build-essential ca-certificates clang cmake curl git jq
  libgmp-dev libmpfr-dev libpq-dev libssl-dev pkg-config ripgrep
  python3 python3-pip python3-venv python3-sympy python3-pytest
  sagemath pari-gp
)
if ! command -v go >/dev/null 2>&1; then
  packages+=(golang-go)
fi

missing=()
for package in "${packages[@]}"; do
  if ! dpkg-query -W -f='${Status}' "$package" 2>/dev/null |
       grep -qx 'install ok installed'; then
    missing+=("$package")
  fi
done
if (("${#missing[@]}" > 0)); then
  run_root env DEBIAN_FRONTEND=noninteractive apt-get update
  run_root env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends "${missing[@]}"
fi

if ! command -v cargo >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -fsSL https://sh.rustup.rs |
    sh -s -- -y --profile minimal
  export PATH="$HOME/.cargo/bin:$PATH"
fi

if ! command -v conductor >/dev/null 2>&1 ||
   ! command -v conductor-mcp >/dev/null 2>&1 ||
   ! command -v conductord >/dev/null 2>&1; then
  GOBIN="$bin_dir" GOFLAGS=-modcacherw     go install "github.com/aburan28/conductor/cmd/...@${CONDUCTOR_VERSION:-main}"
fi

if ! command -v cairn >/dev/null 2>&1; then
  cairn_args=(--locked --root "$HOME/.local" --git https://github.com/aburan28/cairn)
  if [[ -n "${CAIRN_REV:-}" ]]; then
    cairn_args+=(--rev "$CAIRN_REV")
  fi
  cargo install "${cairn_args[@]}" cairn
fi

# Login shells launched after setup can find the user-installed binaries.
profile_line='export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"'
touch "$HOME/.profile"
if ! grep -Fqx "$profile_line" "$HOME/.profile"; then
  printf '%s\n' "$profile_line" >> "$HOME/.profile"
fi

if [[ -f "$repo_dir/Cargo.toml" ]]; then
  (cd "$repo_dir" && cargo fetch --locked)
  # The current main branch has unrelated library compile errors. Enable this
  # after those are fixed to install the repository's Rust executables.
  if [[ "${CRYPTO_BUILD_BINS:-0}" == 1 ]]; then
    cargo install --locked --path "$repo_dir" --bins --root "$HOME/.local"
  fi
fi

echo "Installed tools:"
sage --version
go version
cargo --version
conductor version
cairn --version

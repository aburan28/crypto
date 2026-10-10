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
  build-essential bzip2 ca-certificates clang cmake curl git jq
  libgmp-dev libmpfr-dev libpq-dev libssl-dev pkg-config ripgrep
  python3 python3-pip python3-venv python3-sympy python3-pytest
  pari-gp shellcheck git-lfs
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

# Some minimal images do not expose the APT sagemath package. Prefer it when
# available, then fall back to the Sage project's conda-forge installation.
sage_works() {
  command -v sage >/dev/null 2>&1 &&
    sage -python -c 'from sage.all import GF; assert GF(256).cardinality() == 256' >/dev/null 2>&1
}

if ! sage_works; then
  run_root env DEBIAN_FRONTEND=noninteractive apt-get update
  if apt-cache show sagemath >/dev/null 2>&1; then
    run_root env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends sagemath || true
  fi
fi

if ! sage_works; then
  case "$(uname -m)" in
    x86_64) mamba_platform=linux-64 ;;
    aarch64|arm64) mamba_platform=linux-aarch64 ;;
    *) echo "No micromamba build selected for architecture $(uname -m)" >&2; exit 1 ;;
  esac
  if ! command -v micromamba >/dev/null 2>&1; then
    (
      download_dir="$(mktemp -d)"
      trap 'rm -rf "$download_dir"' EXIT
      curl -fsSL "https://micro.mamba.pm/api/micromamba/$mamba_platform/latest" -o "$download_dir/micromamba.tar.bz2"
      tar -xjf "$download_dir/micromamba.tar.bz2" -C "$download_dir" bin/micromamba
      install -m 755 "$download_dir/bin/micromamba" "$bin_dir/micromamba"
    )
  fi
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.local/share/micromamba}"
  sage_prefix="$MAMBA_ROOT_PREFIX/envs/sage"
  if [[ ! -x "$sage_prefix/bin/sage" ]]; then
    micromamba create -y -p "$sage_prefix" -c conda-forge sage
  fi
  ln -sfn "$sage_prefix/bin/sage" "$bin_dir/sage"
fi

if ! sage_works; then
  echo "SageMath installation finished without a working sage command" >&2
  exit 1
fi

if ! command -v cargo >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -fsSL https://sh.rustup.rs |
    sh -s -- -y --profile minimal
  export PATH="$HOME/.cargo/bin:$PATH"
fi
if command -v rustup >/dev/null 2>&1; then
  rustup component add rustfmt clippy
fi

if ! command -v conductor >/dev/null 2>&1 ||
   ! command -v conductor-mcp >/dev/null 2>&1 ||
   ! command -v conductord >/dev/null 2>&1; then
  GOBIN="$bin_dir" GOFLAGS=-modcacherw go install \
    "github.com/aburan28/conductor/cmd/...@${CONDUCTOR_VERSION:-main}"
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

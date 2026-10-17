#!/usr/bin/env bash
# One-shot workstation setup for the aburan28 research stack:
#
#   * toolchains   - C/C++ (gcc, clang, cmake), Rust (rustup), Go >= 1.25,
#                    Node 22, Python 3 + uv, PostgreSQL 16, PARI/GP, FLINT, GMP
#   * coordination - conductor / conductor-mcp / conductord, cairn
#   * SageMath     - 10.9 from conda-forge via micromamba, plus the
#                    `sage -python` shim that frozen experiments call
#   * agent CLIs   - claude (Claude Code), codex, opencode, gemini (optional)
#   * repositories - crypto, cryptanalysis, crypto-autoresearcher, cloned into
#                    $WORKSPACE and set up (Python venvs, kb, C library build)
#
# Ubuntu/Debian (apt) and macOS (Homebrew). Idempotent: anything already
# present is left alone, so re-running is a cheap way to repair a machine.
# Nothing here starts an autonomous agent (`conductor worker`, `cairn agent run`).
#
#   scripts/setup-workstation.sh                 # everything
#   scripts/setup-workstation.sh --no-sage       # skip the ~8 GiB SageMath env
#   scripts/setup-workstation.sh --only agents   # one component
#   scripts/setup-workstation.sh --list          # component names
#   scripts/setup-workstation.sh --dry-run       # print, change nothing
#
# Environment overrides:
#   WORKSPACE           where repositories are cloned          (~/src)
#   GITHUB_OWNER        owner of the repositories              (aburan28)
#   CONDUCTOR_VERSION   go module version of conductor         (main)
#   CAIRN_REV           cairn git revision                     (default branch)
#   SAGE_VERSION        SageMath version                       (10.9)
#   SAGE_PREFIX         micromamba root for the Sage env       (/opt/conda-sage, or ~/.local/conda-sage without root)
#   GO_VERSION          Go installed when go is missing/old    (1.26.9)
#   NODE_MAJOR          Node.js major version                  (22)
#   WITH_GEMINI=1       also install @google/gemini-cli
#   WITH_LATEX=1        also install a TeX Live (papers/, iacr-eprint skill)
#   WITH_RUST_BUILD=1   also `cargo build --release` the crypto workspace (slow)
set -uo pipefail   # no -e: one failed component warns and the rest still run

case "${WORKSPACE_OFFLINE:-0}" in
  0|'') ;;
  1) exec bash "$(cd "$(dirname "$0")/.." && pwd)/scripts/private-workspace.sh" check ;;
  *) echo "WORKSPACE_OFFLINE must be 0 or 1" >&2; exit 2 ;;
esac

ALL_COMPONENTS=(system rust go node python conductor cairn sage agents repos control-plane)

WORKSPACE="${WORKSPACE:-$HOME/src}"
GITHUB_OWNER="${GITHUB_OWNER:-aburan28}"
CONDUCTOR_VERSION="${CONDUCTOR_VERSION:-main}"
SAGE_VERSION="${SAGE_VERSION:-10.9}"
GO_VERSION="${GO_VERSION:-1.26.9}"
GO_MIN="1.25"
NODE_MAJOR="${NODE_MAJOR:-22}"
REPOS=(crypto cryptanalysis crypto-autoresearcher)

DRY_RUN=0
declare -A SKIP=()
ONLY=()

usage() { sed -n '2,35p' "$0" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY_RUN=1 ;;
    --only) shift; IFS=, read -r -a ONLY <<<"${1:?--only needs a component list}" ;;
    --skip) shift; IFS=, read -r -a s <<<"${1:?--skip needs a component list}"; for c in "${s[@]}"; do SKIP[$c]=1; done ;;
    --no-sage) SKIP[sage]=1 ;;
    --no-agents) SKIP[agents]=1 ;;
    --no-repos) SKIP[repos]=1 ;;
    --no-control-plane) SKIP[control-plane]=1 ;;
    --workspace) shift; WORKSPACE="${1:?}" ;;
    --list) printf '%s\n' "${ALL_COMPONENTS[@]}"; exit 0 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
  esac
  shift
done

# ---------------------------------------------------------------- helpers ---
log()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
ok()   { printf '\033[1;32m ok\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m !!\033[0m %s\n' "$*" >&2; WARNINGS+=("$*"); }
WARNINGS=()
have() { command -v "$1" >/dev/null 2>&1; }
run()  { if [ "$DRY_RUN" = 1 ]; then printf '   + %s\n' "$*"; else "$@"; fi; }
want() {
  local c="$1"
  [ -n "${SKIP[$c]:-}" ] && return 1
  [ "${#ONLY[@]}" -eq 0 ] && return 0
  local o; for o in "${ONLY[@]}"; do [ "$o" = "$c" ] && return 0; done
  return 1
}
# version_ge A B: true when A >= B (dotted numeric).
version_ge() { [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]; }

OS="$(uname -s)"
ARCH="$(uname -m)"
case "$ARCH" in x86_64|amd64) GOARCH=amd64; MMARCH=64 ;; aarch64|arm64) GOARCH=arm64; MMARCH=aarch64 ;; *) GOARCH="$ARCH"; MMARCH="$ARCH" ;; esac

if [ "$(id -u)" -eq 0 ]; then SUDO=(); elif have sudo; then SUDO=(sudo); else SUDO=(); fi
as_root() { if [ "$(id -u)" -eq 0 ] || [ "${#SUDO[@]}" -gt 0 ]; then run "${SUDO[@]}" "$@"; else warn "needs root, skipped: $*"; return 1; fi; }

PKG=""
if [ "$OS" = Darwin ]; then PKG=brew; elif have apt-get; then PKG=apt; elif have dnf; then PKG=dnf; fi

LOCAL_BIN="$HOME/.local/bin"
mkdir -p "$LOCAL_BIN"
export PATH="$LOCAL_BIN:$HOME/.cargo/bin:/usr/local/go/bin:$HOME/go/bin:$PATH"

# Every PATH entry this script relies on, written once to the user's shell rc.
persist_path() {
  local line="$1" rc
  for rc in "$HOME/.bashrc" "$HOME/.zshrc"; do
    [ -f "$rc" ] || [ "$rc" = "$HOME/.bashrc" ] || continue
    grep -qxF "$line" "$rc" 2>/dev/null || run sh -c "printf '%s\n' '$line' >> '$rc'"
  done
}

# ----------------------------------------------------------------- system ---
install_system() {
  log "system packages ($PKG)"
  case "$PKG" in
    apt)
      local pkgs=(build-essential clang lld cmake ninja-build pkg-config git curl wget ca-certificates
        unzip xz-utils bzip2 jq ripgrep shellcheck cppcheck valgrind clang-tidy
        python3 python3-dev python3-pip python3-venv
        libgmp-dev libmpfr-dev libmpc-dev libflint-dev libssl-dev libffi-dev
        pari-gp libpari-dev postgresql postgresql-contrib libpq-dev
        poppler-utils ghostscript graphviz)
      [ "${WITH_LATEX:-0}" = 1 ] && pkgs+=(texlive-latex-extra texlive-science texlive-fonts-recommended latexmk)
      as_root env DEBIAN_FRONTEND=noninteractive apt-get update -qq || true
      as_root env DEBIAN_FRONTEND=noninteractive apt-get install -y -qq --no-install-recommends "${pkgs[@]}" \
        || warn "some apt packages failed; re-run or install them by hand"
      ;;
    dnf)
      as_root dnf install -y gcc gcc-c++ clang cmake ninja-build git curl jq ripgrep ShellCheck cppcheck valgrind \
        python3 python3-devel gmp-devel mpfr-devel libmpc-devel flint-devel openssl-devel pari-gp \
        postgresql-server postgresql-contrib libpq-devel poppler-utils graphviz \
        || warn "some dnf packages failed"
      ;;
    brew)
      have brew || { warn "Homebrew not found: install it from https://brew.sh and re-run"; return; }
      run brew install cmake ninja pkg-config git jq ripgrep shellcheck cppcheck llvm gmp mpfr libmpc flint \
        pari postgresql@16 poppler graphviz python@3.12 || warn "some brew packages failed"
      [ "${WITH_LATEX:-0}" = 1 ] && run brew install --cask mactex-no-gui
      ;;
    *) warn "no supported package manager; install a C toolchain, cmake, GMP/MPFR/FLINT, PARI/GP and PostgreSQL 16 yourself" ;;
  esac
}

# ------------------------------------------------------------------- rust ---
install_rust() {
  log "Rust"
  if ! have rustup; then
    run sh -c 'curl --proto =https --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile default'
    export PATH="$HOME/.cargo/bin:$PATH"
  fi
  run rustup toolchain install stable --profile default --component clippy,rustfmt
  run rustup default stable
  have cargo-deny || run cargo install --locked cargo-deny || warn "cargo-deny not installed (cryptanalysis 'make suite-lint')"
  persist_path 'export PATH="$HOME/.cargo/bin:$PATH"'
  have cargo && ok "$(cargo --version)"
}

# --------------------------------------------------------------------- go ---
install_go() {
  log "Go (>= $GO_MIN)"
  local cur=""
  have go && cur="$(go env GOVERSION 2>/dev/null | sed 's/^go//')"
  if [ -n "$cur" ] && version_ge "$cur" "$GO_MIN"; then ok "go $cur"; return; fi
  if [ "$PKG" = brew ]; then run brew install go; return; fi
  local os=linux; [ "$OS" = Darwin ] && os=darwin
  local tgz="go${GO_VERSION}.${os}-${GOARCH}.tar.gz"
  run sh -c "curl -fsSL https://go.dev/dl/$tgz -o /tmp/$tgz"
  as_root rm -rf /usr/local/go && as_root tar -C /usr/local -xzf "/tmp/$tgz" && run rm -f "/tmp/$tgz"
  persist_path 'export PATH="/usr/local/go/bin:$HOME/go/bin:$PATH"'
  have go && ok "$(go version)"
}

# ------------------------------------------------------------------- node ---
install_node() {
  log "Node.js $NODE_MAJOR"
  if have node && [ "$(node -p 'process.versions.node.split(".")[0]')" -ge "$NODE_MAJOR" ]; then ok "node $(node --version)"; return; fi
  case "$PKG" in
    brew) run brew install "node@$NODE_MAJOR" ;;
    apt) run sh -c "curl -fsSL https://deb.nodesource.com/setup_${NODE_MAJOR}.x | ${SUDO[*]:-} bash -" \
           && as_root apt-get install -y -qq nodejs ;;
    *) warn "install Node.js >= $NODE_MAJOR yourself (https://nodejs.org)" ;;
  esac
}

# Global npm installs into ~/.local so they never need sudo.
npm_global() {
  have npm || { warn "npm missing, cannot install $*"; return 1; }
  if [ ! -w "$(npm prefix -g 2>/dev/null)/lib" ] 2>/dev/null; then
    run npm config set prefix "$HOME/.local"
  fi
  run npm install -g "$@"
}

# ----------------------------------------------------------------- python ---
install_python() {
  log "Python tooling (uv + research libraries)"
  have uv || run sh -c 'curl -LsSf https://astral.sh/uv/install.sh | sh'
  export PATH="$HOME/.local/bin:$PATH"
  persist_path 'export PATH="$HOME/.local/bin:$PATH"'
  # User-level libraries the plotting/report scripts in crypto/ import
  # (matplotlib, reportlab, numpy, ...). Repo-specific deps go in each repo's venv.
  run uv tool install --quiet ruff || true
  run uv tool install --quiet pre-commit || true
  have uv && ok "$(uv --version)"
}

# -------------------------------------------------------------- conductor ---
install_conductor() {
  log "Conductor ($CONDUCTOR_VERSION)"
  if have conductor && have conductor-mcp && have conductord; then ok "$(conductor version 2>/dev/null | head -n1)"; return; fi
  have go || { warn "go missing; run with --only go first"; return 1; }
  # @latest resolves to v0.1.0 whose go.mod has the wrong module path: use @main.
  run env GOBIN="$LOCAL_BIN" GOFLAGS=-modcacherw \
    go install "github.com/aburan28/conductor/cmd/...@${CONDUCTOR_VERSION}" \
    || warn "conductor install failed (try: curl -fsSL https://raw.githubusercontent.com/aburan28/conductor/main/scripts/install-release.sh | bash -s -- aburan28/conductor ~/.local/bin)"
}

# ------------------------------------------------------------------ cairn ---
install_cairn() {
  log "cairn"
  if have cairn; then ok "$(cairn --version 2>/dev/null | head -n1)"; return; fi
  have cargo || { warn "cargo missing; run with --only rust first"; return 1; }
  run cargo install --locked --root "$HOME/.local" \
    --git "https://github.com/$GITHUB_OWNER/cairn" ${CAIRN_REV:+--rev "$CAIRN_REV"} cairn \
    || warn "cairn install failed"
  mkdir -p "$HOME/.cairn"
}

# ------------------------------------------------------------------- sage ---
install_sage() {
  log "SageMath $SAGE_VERSION (conda-forge, ~8 GiB)"
  if [ -z "${SAGE_PREFIX:-}" ]; then
    if [ -w /opt ] || [ "${#SUDO[@]}" -gt 0 ]; then SAGE_PREFIX=/opt/conda-sage; else SAGE_PREFIX="$HOME/.local/conda-sage"; fi
  fi
  local env="$SAGE_PREFIX/envs/sage"
  if [ -x "$env/bin/python3" ] && "$env/bin/python3" -c \
      "from sage.env import SAGE_VERSION as v; import sys; sys.exit(v != '$SAGE_VERSION')" 2>/dev/null; then
    ok "sage $SAGE_VERSION at $env"
  else
    local parent; parent="$(dirname "$SAGE_PREFIX")"
    local free_kb; free_kb="$(df -Pk "$parent" 2>/dev/null | awk 'NR==2{print $4}')"
    if [ -n "$free_kb" ] && [ "$free_kb" -lt 10485760 ]; then
      warn "<10 GiB free under $parent, skipping SageMath (set SAGE_PREFIX elsewhere)"; return
    fi
    if [ ! -d "$SAGE_PREFIX" ]; then
      if [ -w "$parent" ]; then run mkdir -p "$SAGE_PREFIX"; else as_root install -d -o "$(id -u)" -g "$(id -g)" "$SAGE_PREFIX"; fi
    fi
    local mm="$LOCAL_BIN/micromamba"
    if ! have micromamba && [ ! -x "$mm" ]; then
      local plat=linux; [ "$OS" = Darwin ] && plat=osx
      run sh -c "curl -fsSL https://micro.mamba.pm/api/micromamba/${plat}-${MMARCH}/latest | tar xj -C '$HOME/.local' bin/micromamba"
    fi
    have micromamba && mm="$(command -v micromamba)"
    run "$mm" create -y -r "$SAGE_PREFIX" -p "$env" -c conda-forge "sage=$SAGE_VERSION" \
      || { warn "SageMath install failed"; return; }
  fi
  # `sage -python script.py` is the upstream CLI convention frozen trial plans
  # use; conda-forge's `sage` has no -python, so wrap it.
  run mkdir -p "$SAGE_PREFIX/bin"
  if [ "$DRY_RUN" = 0 ]; then
    cat >"$SAGE_PREFIX/bin/sage" <<EOF
#!/bin/sh
if [ "\$1" = "-python" ]; then shift; exec "$env/bin/python3" "\$@"; fi
exec "$env/bin/sage" "\$@"
EOF
    chmod +x "$SAGE_PREFIX/bin/sage"
  fi
  persist_path "export PATH=\"$SAGE_PREFIX/bin:\$PATH\""
  export PATH="$SAGE_PREFIX/bin:$PATH"
}

# ----------------------------------------------------------------- agents ---
install_agents() {
  log "agent CLIs"
  if have claude; then ok "claude $(claude --version 2>/dev/null | head -n1)"
  else run sh -c 'curl -fsSL https://claude.ai/install.sh | bash' || npm_global @anthropic-ai/claude-code || warn "claude not installed"; fi
  if have codex; then ok "codex $(codex --version 2>/dev/null | head -n1)"
  else npm_global @openai/codex || warn "codex not installed"; fi
  if have opencode; then ok "opencode $(opencode --version 2>/dev/null | head -n1)"
  else npm_global opencode-ai || run sh -c 'curl -fsSL https://opencode.ai/install | bash' || warn "opencode not installed"; fi
  if [ "${WITH_GEMINI:-0}" = 1 ]; then have gemini || npm_global @google/gemini-cli || warn "gemini not installed"; fi
  have gh || case "$PKG" in brew) run brew install gh ;; apt) as_root apt-get install -y -qq gh || warn "gh not installed" ;; esac
}

# ------------------------------------------------------------------ repos ---
clone_repo() {
  local name="$1" dir="$WORKSPACE/$1"
  if [ -d "$dir/.git" ]; then ok "$name present at $dir"; return; fi
  run git clone "https://github.com/$GITHUB_OWNER/$name.git" "$dir" || { warn "clone of $name failed"; return 1; }
}

setup_crypto() {
  local dir="$WORKSPACE/crypto"; [ -d "$dir" ] || return 0
  log "crypto: Python venv for plotting/report scripts"
  run uv venv --quiet --allow-existing "$dir/.venv"
  run uv pip install --quiet --python "$dir/.venv/bin/python" \
    numpy scipy sympy matplotlib pandas reportlab pyyaml networkx cypari2 python-flint pytest \
    || warn "crypto venv: some packages failed"
  [ "${WITH_RUST_BUILD:-0}" = 1 ] && { log "crypto: cargo build --release"; (cd "$dir" && run cargo build --release) || warn "crypto build failed"; }
  return 0
}

setup_cryptanalysis() {
  local dir="$WORKSPACE/cryptanalysis"; [ -d "$dir" ] || return 0
  log "cryptanalysis: C library (cmake) + Python bindings"
  (cd "$dir" && run make lib) || warn "cryptanalysis: make lib failed"
  run uv venv --quiet --allow-existing "$dir/.venv"
  run uv pip install --quiet --python "$dir/.venv/bin/python" ruff pytest numpy || true
  [ -f "$dir/bindings/python/pyproject.toml" ] && \
    { run uv pip install --quiet --python "$dir/.venv/bin/python" -e "$dir/bindings/python" || warn "cryptanalysis python bindings failed"; }
  [ -f "$dir/suite/python/indexcalc/requirements-sat.txt" ] && \
    { run uv pip install --quiet --python "$dir/.venv/bin/python" -r "$dir/suite/python/indexcalc/requirements-sat.txt" || true; }
  return 0
}

setup_autoresearcher() {
  local dir="$WORKSPACE/crypto-autoresearcher"; [ -d "$dir" ] || return 0
  log "crypto-autoresearcher: editable install + kb"
  run uv venv --quiet --allow-existing "$dir/.venv"
  run uv pip install --quiet --python "$dir/.venv/bin/python" \
    -r "$dir/requirements-dev.txt" -r "$dir/requirements-agent.txt" \
    -e "${dir}[dev,fast,gf2,index-calculus,agent,research-loop,campaign-mcp,rho-dp]" \
    || warn "crypto-autoresearcher: install failed"
  [ -d "$dir/kb" ] && { (cd "$dir/kb" && run uv sync --quiet --extra index --extra mcp --extra dev) || warn "kb: uv sync failed"; }
  [ -d "$dir/.git" ] && { (cd "$dir" && run make hooks) || true; }
  return 0
}

install_repos() {
  log "repositories in $WORKSPACE"
  run mkdir -p "$WORKSPACE"
  local r; for r in "${REPOS[@]}"; do clone_repo "$r" || true; done
  setup_crypto; setup_cryptanalysis; setup_autoresearcher
}

# ---------------------------------------------------------- control plane ---
start_control_plane() {
  log "Conductor control plane"
  have conductor || { warn "conductor missing, no control plane"; return; }
  if [ -n "${CONDUCTOR_ENDPOINT:-}" ] && [ -n "${CONDUCTOR_TOKEN:-}" ]; then
    ok "using shared control plane at $CONDUCTOR_ENDPOINT"; return
  fi
  local dir="$WORKSPACE/crypto"; [ -d "$dir" ] || dir="$PWD"
  # `conductor up` uses a reachable Postgres, else Docker, and prints a login token.
  (cd "$dir" && run conductor up --project crypto) || warn "conductor up failed; see docs/CONDUCTOR.md"
  (cd "$dir" && run conductor doctor) || true
}

# ------------------------------------------------------------------- main ---
log "workstation setup: $OS/$ARCH, package manager: ${PKG:-none}, workspace: $WORKSPACE"
want system        && install_system || true
want rust          && install_rust || true
want go            && install_go || true
want node          && install_node || true
want python        && install_python || true
want conductor     && install_conductor || true
want cairn         && install_cairn || true
want sage          && install_sage || true
want agents        && install_agents || true
want repos         && install_repos || true
want control-plane && start_control_plane || true

echo
log "summary"
for c in gcc clang cmake cargo go node uv psql gp conductor conductor-mcp cairn sage claude codex opencode gh; do
  if have "$c"; then printf '   %-14s %s\n' "$c" "$(command -v "$c")"; else printf '   %-14s \033[1;33mmissing\033[0m\n' "$c"; fi
done
if [ "${#WARNINGS[@]}" -gt 0 ]; then
  echo; log "warnings"; printf '   - %s\n' "${WARNINGS[@]}"
fi
echo
echo "Open a new shell (or 'source ~/.bashrc') to pick up PATH changes."
echo "API keys are not handled here: set ANTHROPIC_API_KEY / OPENAI_API_KEY or run 'claude', 'codex login', 'opencode auth login'."

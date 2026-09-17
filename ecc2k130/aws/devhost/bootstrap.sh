#!/usr/bin/env bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
export NEEDRESTART_MODE=a

log() { printf '[g7e-devhost] %s\n' "$*"; }

log "installing base packages"
apt-get update -y
apt-get install -y ca-certificates curl git jq build-essential tmux htop unzip ripgrep fd-find python3 python3-pip python3-venv

# Node 22 is used for the coding CLIs. Re-running this script upgrades in place.
if ! command -v node >/dev/null 2>&1 || [ "$(node -p 'Number(process.versions.node.split(`.`)[0])' 2>/dev/null || echo 0)" -lt 22 ]; then
  curl -fsSL https://deb.nodesource.com/setup_22.x | bash -
  apt-get install -y nodejs
fi

log "installing GitHub CLI"
install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg -o /etc/apt/keyrings/githubcli-archive-keyring.gpg
chmod go+r /etc/apt/keyrings/githubcli-archive-keyring.gpg
printf 'deb [arch=%s signed-by=/etc/apt/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main\n' "$(dpkg --print-architecture)" > /etc/apt/sources.list.d/github-cli.list
apt-get update -y
apt-get install -y gh

log "installing/upgrading Claude Code, Codex, and OpenCode"
npm install -g @anthropic-ai/claude-code@latest @openai/codex@latest opencode-ai@latest

# Durable working area lives on EBS. Instance-store NVMe is intentionally not used
# for anything that must survive stop/start.
mkdir -p /workspace
chmod 0777 /workspace

cat >/usr/local/bin/refresh-agent-tools <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
npm install -g @anthropic-ai/claude-code@latest @openai/codex@latest opencode-ai@latest
apt-get update -y >/dev/null
apt-get install -y gh >/dev/null
printf 'claude:  '; claude --version || true
printf 'codex:   '; codex --version || true
printf 'opencode:'; opencode --version || true
printf 'gh:      '; gh --version | head -1 || true
EOF
chmod 0755 /usr/local/bin/refresh-agent-tools

cat >/etc/systemd/system/g7e-agent-tools-refresh.service <<'EOF'
[Unit]
Description=Refresh coding-agent CLIs after boot
After=network-online.target
Wants=network-online.target

[Service]
Type=oneshot
ExecStart=/usr/local/bin/refresh-agent-tools

[Install]
WantedBy=multi-user.target
EOF
systemctl daemon-reload
systemctl enable g7e-agent-tools-refresh.service

# Helpful shell defaults for interactive SSH/SSM sessions.
cat >/etc/profile.d/g7e-devhost.sh <<'EOF'
export WORKSPACE=/workspace
export PATH=/usr/local/bin:$PATH
EOF
chmod 0644 /etc/profile.d/g7e-devhost.sh

log "tool versions"
node --version
npm --version
gh --version | head -1 || true
claude --version || true
codex --version || true
opencode --version || true
nvidia-smi || true

# Runaway cost controls: optional host timer (idle/age default OFF).
# Monthly $5k enforcement lives in costguard/costguard.sh on a control host.
log "installing costguard host timer (idle/age disabled by default)"
install -d -m 755 /var/lib/costguard /usr/local/lib/costguard
if [[ -x /home/ubuntu/crypto/ecc2k130/aws/costguard/install-host.sh ]]; then
  COSTGUARD_MAX_AGE_HOURS=${COSTGUARD_MAX_AGE_HOURS:-0} \
  COSTGUARD_IDLE_HOURS=${COSTGUARD_IDLE_HOURS:-0} \
    bash /home/ubuntu/crypto/ecc2k130/aws/costguard/install-host.sh
elif [[ -x /workspace/crypto/ecc2k130/aws/costguard/install-host.sh ]]; then
  COSTGUARD_MAX_AGE_HOURS=${COSTGUARD_MAX_AGE_HOURS:-0} \
  COSTGUARD_IDLE_HOURS=${COSTGUARD_IDLE_HOURS:-0} \
    bash /workspace/crypto/ecc2k130/aws/costguard/install-host.sh
else
  src=/tmp/costguard-src
  mkdir -p "$src" \
    && curl -fsSL "https://raw.githubusercontent.com/aburan28/crypto/main/ecc2k130/aws/costguard/install-host.sh" \
      -o "$src/install-host.sh" \
    && curl -fsSL "https://raw.githubusercontent.com/aburan28/crypto/main/ecc2k130/aws/costguard/host-idle-stop.sh" \
      -o "$src/host-idle-stop.sh" \
    && COSTGUARD_MAX_AGE_HOURS=0 COSTGUARD_IDLE_HOURS=0 bash "$src/install-host.sh" \
    || log "costguard install deferred (scripts not yet on main / no checkout)"
fi

log "bootstrap complete; authenticate interactively with: gh auth login, claude, codex, and opencode"

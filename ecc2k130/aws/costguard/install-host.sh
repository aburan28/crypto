#!/usr/bin/env bash
#
# Install optional host-local costguard timer.
# Defaults: idle/age DISABLED (monthly budget is the primary control).
#
#   sudo ./install-host.sh
#   COSTGUARD_MAX_AGE_HOURS=168 COSTGUARD_IDLE_HOURS=48 sudo -E ./install-host.sh
#
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MAX_AGE_HOURS=${COSTGUARD_MAX_AGE_HOURS:-0}
IDLE_HOURS=${COSTGUARD_IDLE_HOURS:-0}
IDLE_LOAD=${COSTGUARD_IDLE_LOAD:-0.3}

install -d -m 755 /var/lib/costguard /usr/local/lib/costguard
install -m 755 "$ROOT/host-idle-stop.sh" /usr/local/lib/costguard/host-idle-stop.sh

cat > /etc/default/costguard <<EOF
COSTGUARD_MAX_AGE_HOURS=$MAX_AGE_HOURS
COSTGUARD_IDLE_HOURS=$IDLE_HOURS
COSTGUARD_IDLE_LOAD=$IDLE_LOAD
COSTGUARD_READYZ_URL=http://127.0.0.1:9182/readyz
EOF
chmod 644 /etc/default/costguard

cat > /etc/systemd/system/costguard-idle-stop.service <<'EOF'
[Unit]
Description=Optional GPU host idle/age cost check
After=network-online.target

[Service]
Type=oneshot
EnvironmentFile=-/etc/default/costguard
ExecStart=/usr/local/lib/costguard/host-idle-stop.sh
EOF

cat > /etc/systemd/system/costguard-idle-stop.timer <<'EOF'
[Unit]
Description=Run optional costguard idle/age check every 15 minutes

[Timer]
OnBootSec=10min
OnUnitActiveSec=15min
AccuracySec=1min
Unit=costguard-idle-stop.service

[Install]
WantedBy=timers.target
EOF

date +%s > /var/lib/costguard/armed-at
touch /var/lib/costguard/keepalive
chmod 644 /var/lib/costguard/armed-at /var/lib/costguard/keepalive

systemctl daemon-reload
systemctl enable --now costguard-idle-stop.timer
echo "costguard host timer enabled (max_age=${MAX_AGE_HOURS}h idle=${IDLE_HOURS}h; 0=disabled)"
systemctl list-timers costguard-idle-stop.timer --no-pager || true

#!/usr/bin/env bash
# Install dp_ingest.py as a per-user LaunchAgent.
#
# Copies the ingest scripts into ~/Library/Application Support and bootstraps
# com.adamburan.ecc2k130-dp-ingest. The agent starts at login. launchd restarts
# it when the process exits non-zero. Each start refreshes this machine's
# egress /32 on the rho-dp security group (INGEST_ENSURE_ACCESS=1).
#
# The running agent uses the copies made here. After a new dp_ingest.py lands,
# run this script again.
set -euo pipefail

launchd_dir=$(cd "$(dirname "$0")" && pwd)
aws_dir=$(cd "$launchd_dir/.." && pwd)
home=${HOME:?}
label=com.adamburan.ecc2k130-dp-ingest
ingest="$home/Library/Application Support/ECC2K130/ingest"
logs="$home/Library/Logs/ECC2K130"
agents="$home/Library/LaunchAgents"
plist="$agents/$label.plist"
uid=$(id -u)
domain="gui/$uid"

mkdir -p "$ingest" "$logs" "$agents"
cp "$aws_dir/dp_ingest.py" "$aws_dir/ingest.sh" "$launchd_dir/run-dp-ingest.zsh" "$ingest/"
chmod 755 "$ingest/ingest.sh" "$ingest/dp_ingest.py" "$ingest/run-dp-ingest.zsh"

cat > "$plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key>
  <string>$label</string>

  <key>ProgramArguments</key>
  <array>
    <string>/bin/zsh</string>
    <string>$ingest/run-dp-ingest.zsh</string>
  </array>

  <key>WorkingDirectory</key>
  <string>$ingest</string>

  <key>RunAtLoad</key>
  <true/>

  <key>KeepAlive</key>
  <dict>
    <key>SuccessfulExit</key>
    <false/>
  </dict>

  <key>ExitTimeOut</key>
  <integer>120</integer>

  <key>ProcessType</key>
  <string>Background</string>

  <key>ThrottleInterval</key>
  <integer>30</integer>

  <key>StandardOutPath</key>
  <string>$logs/dp-ingest.stdout.log</string>

  <key>StandardErrorPath</key>
  <string>$logs/dp-ingest.stderr.log</string>
</dict>
</plist>
EOF

launchctl bootout "$domain/$label" 2>/dev/null || true
launchctl bootstrap "$domain" "$plist"
launchctl enable "$domain/$label"
launchctl kickstart -k "$domain/$label"
echo "bootstrapped $domain/$label"

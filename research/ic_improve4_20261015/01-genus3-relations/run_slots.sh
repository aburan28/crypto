#!/bin/sh
# Keep at most $SLOTS (default 3) of this item's jobs running, taking commands (job.sh
# arguments) from a queue file in order. Usage: sh run_slots.sh queue.txt
Q=$1
running() { pgrep -fc '^\.\./\.\./\.\./target/release/examples/ic_improve4_01_genus3_rel' ; }
while read -r line; do
  case "$line" in ''|'#'*) continue ;; esac
  while [ "$(running)" -ge "${SLOTS:-3}" ]; do sleep 10; done
  sh job.sh $line &
  sleep 3
done < "$Q"
wait

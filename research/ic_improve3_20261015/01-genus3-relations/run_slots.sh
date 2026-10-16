#!/bin/sh
# Keep at most 2 of this item's jobs running (counting jobs started by any script), taking
# commands (job.sh arguments) from a queue file in order. Usage: sh run_slots.sh queue.txt
Q=$1
running() { pgrep -fc '^\.\./\.\./\.\./target/release/examples/ic_improve3_01_genus3_rel' ; }
while read -r line; do
  while [ "$(running)" -ge 2 ]; do sleep 10; done
  sh job.sh $line &
  sleep 5
done < "$Q"
wait

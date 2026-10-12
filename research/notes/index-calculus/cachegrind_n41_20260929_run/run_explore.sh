#!/bin/bash
# Amendment A1 follow-ups, quiet machine, in this order.
set -u
cd "$(dirname "$0")"
: "${MEM_CALIB:?set MEM_CALIB}"
for st in replicate thp interference; do
  echo "== explore $st start $(date -u +%FT%TZ) load: $(cat /proc/loadavg)" >> stages.log
  python3 explore_stages.py $st > work/explore_$st.log 2>&1
  echo "== explore $st rc=$? end $(date -u +%FT%TZ)" >> stages.log
done
echo EXPLORE_DONE >> stages.log

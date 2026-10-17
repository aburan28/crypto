#!/bin/sh
# wait for the diagnosis queue, then the test priorities 1-2 (sched, curves 0-2, 12..21; then b = 22)
while pgrep -f "^sh run_diag.sh" >/dev/null; do sleep 20; done
PAR=2 sh run_test.sh results/run sched "12 13 14 15 16 17 18 19 20 21" "0 1 2"
PAR=2 sh run_test.sh results/run sched "22" "0 1 2"

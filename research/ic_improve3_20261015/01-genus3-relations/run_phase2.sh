#!/bin/sh
# after phase 1 (sched curves 0-2 incl. b = 22): priorities 3-5 (needs r2.env = amendment C1)
while pgrep -f "^sh run_phase1.sh" >/dev/null; do sleep 20; done
PAR=2 sh run_test.sh results/run r2 "12 13 14 15 16 17 18 19 20 21" "0 1 2"
PAR=2 sh run_test.sh results/run aa "12 13 14 15 16 17 18 19 20 21" "0 1 2"
PAR=2 sh run_test.sh results/run r2 "22" "0 1 2"
PAR=2 sh run_test.sh results/run sched "18 19 20 21" "3 4"
PAR=2 sh run_test.sh results/run r2 "18 19 20 21" "3 4"

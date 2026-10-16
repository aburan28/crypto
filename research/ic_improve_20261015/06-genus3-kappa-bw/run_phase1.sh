#!/bin/sh
# Phase 1 orchestration: training b = 20, 21 (PAR 2) alongside the main/aa test collections (PAR 1).
while pgrep -f "run_train.sh results/train 16 18" >/dev/null; do sleep 5; done
PAR=2 sh run_train.sh results/train "20 21" &
mkdir -p results/run
PAR=1 sh run_test.sh results/run main "12 13 14 15 16 17 18 19 20" "0 1 2"
PAR=1 sh run_test.sh results/run aa "12 13 14 15 16 17 18 19 20" "0 1 2"
wait

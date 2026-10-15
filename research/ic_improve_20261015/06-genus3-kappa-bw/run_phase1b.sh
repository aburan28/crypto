#!/bin/sh
# After the phase-1 main/aa queue (b <= 20): main at b = 21, 22 on curves 0-2 (PAR 1).
until [ -s results/run/b20_c2_aa.jsonl ]; do sleep 10; done
PAR=1 sh run_test.sh results/run main "21 22" "0 1 2"

#!/bin/sh
PAR=2 sh run_test.sh results/run sched "12 13 14 15 16 17 18 19 20 21" "0 1 2"
PAR=2 sh run_test.sh results/run sched "18 19 20 21" "3 4"

#!/bin/sh
# Part D, E and F claim sessions (PROTOCOL.md), four streams on four CPUs,
# each with its own lock.  Counts are deterministic; wall time is
# descriptive (L1: the streams share the host).
cd "$(dirname "$0")/../.."
R=research/gaudry_k4_char2_crossover_20261015
E=target2/release/ecbench
run() { # cpu lock spec out
  $E run --spec $R/specs/$3.json --out $R/sessions/$4 --cpus $1 --lock /tmp/bench-$2.lock --wait --allow-busy > $R/sessions/$4.log 2>&1
}
mkdir -p $R/sessions
( run 1 s1 SPEC-E-m18-rho e-m18-rho ) &
( run 2 s2 SPEC-E-m18-ic e-m18-ic; run 2 s2 SPEC-F-lp f-lp ) &
( run 3 s3 SPEC-E-m17-ic e-m17-ic ) &
( run 0 s0 SPEC-E-m17-rho e-m17-rho; run 0 s0 SPEC-E-m16 e-m16; run 0 s0 SPEC-E-m15 e-m15; run 0 s0 SPEC-D-ladder d-ladder; run 0 s0 SPEC-D-m14 d-m14 ) &
wait

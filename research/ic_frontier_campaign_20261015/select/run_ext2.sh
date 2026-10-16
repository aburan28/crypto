#!/bin/sh
while pgrep -f run_select.sh >/dev/null; do sleep 30; done
for n in 47 41; do
  ../../../target/release/ecbench run --spec spec-select-K$n-ext.json --out sessions/K$n-ext --cpus none --allow-busy --quiet --lock /tmp/ic-frontier-ext.lock
done

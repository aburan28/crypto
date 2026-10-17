# Amendment 3 (2026-10-14, before the re-run data)

Builds of the instrument and diagnostic runs overlapped 38 recorded runs
(`raw/contended_runs.jsonl`, computed from each record's end time and
elapsed time against the windows in `raw/contention.log`, widened by 30 s
where the log has no end time). Each of those (pass, arm, D, m, seed,
threads) cells is re-run once on the idle host with the same binary its
pass used, into `raw/rerun_contended.jsonl`. The report uses the re-run
value for those cells and lists both. No other cell changes.

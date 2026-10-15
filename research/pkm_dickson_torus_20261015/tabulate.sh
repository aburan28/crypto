#!/bin/bash
# Lay out the recorded rows of runs/*.jsonl as one TSV (no computation beyond
# field extraction): cells.tsv, one line per measured system.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
{
  printf 'run\tkind\tp\tm\tN\tcontrol\ttarget_index\tD\twidth_to_solution\tsteps_to_solution\tinconsistent\ttimed_out\toversize\tmuladds\tms\n'
  for f in runs/R*.jsonl; do
    r=$(basename "$f" .jsonl)
    jq -r --arg r "$r" 'select(.kind) | [$r, .kind, .p, .m, .N, .control, .target_index, .solving_degree_max, .max_cols_to_solution, .steps_to_solution, .inconsistent, .timed_out, (.oversize // false), (.muladds // ""), .ms] | @tsv' "$f"
  done
} > cells.tsv
cat > recorded_split.tsv <<'TSV'
source	kind	p	m	N	D	width_to_solution	steps_to_solution	target
round2/XV1	kummer=dickson	786433	2	2	3	4	1	planted
round2/XV1	kummer=dickson	786433	2	4	4	16	3	planted
round2/XV1	kummer=dickson	786433	2	6	4	56	5	planted
round2/XV1	kummer=dickson	786433	2	8	4	152	7	planted
round2/XV1	kummer=dickson	786433	2	10	4	337/336	10	planted
round2/XV1	kummer=dickson	786433	2	12	5	1231/1223	13	planted
round2/XV1	kummer=dickson	786433	2	14	5	2304/2303	15	planted
round2/XV1	kummer=dickson	786433	2	16	5	4789/4786	17	random
round2/XV1	kummer=dickson	786433	2	18	5	8644/8639	20	random
round2/K0	kummer	786433	2	20	6	33446	24	random
round2/K1	kummer	2013265921	2	12	5	1231	13	random
round2/K1	kummer	2013265921	2	14	5	2304	15	random
round2/K1	kummer	2013265921	2	16	5	4789	17	random
round2/K1	kummer	2013265921	2	18	5	8644	20	random
round2/K1	kummer	2013265921	2	20	6	33446	24	random
round2/XV5	kummer	3221225473	2	2	3	4	1	planted
round2/XV5	kummer	3221225473	2	4	4	16	3	planted
round2/XV5	kummer	3221225473	2	6	4	56	5	planted
round2/XV5	kummer	3221225473	2	8	4	152	7	planted
round2/XV5	kummer	3221225473	2	10	4	337	10	planted
round2/XV5	kummer	3221225473	2	12	5	1231	13	planted
round2/XV5	kummer	3221225473	2	14	5	2304	15	planted
round2/XV5	kummer	3221225473	2	16	5	4789	17	random
round2/XV2	kummer=dickson	786433	3	6	5	164	10	planted
round2/XV2	kummer=dickson	786433	3	9	6	1186/1184	16	planted
round2/XV2	kummer=dickson	786433	3	12	6	5019/5018	18	random
round2/M3	kummer	2013265921	3	9	6	1186	16	random
round2/M3	kummer	2013265921	3	12	6	5019	18	random
round2/M3	kummer	2013265921	3	15	7	31397	24	random
round2/XV3	kummer	786433	4	8	6	1049	21	random
round2/XV3	kummer	786433	4	12	7	9452	29	random
round2/M4	kummer	2013265921	4	8	6	1049	21	random
round2/M4	kummer	2013265921	4	12	7	9452	29	random
TSV
echo "wrote cells.tsv ($(($(wc -l < cells.tsv) - 1)) systems) and recorded_split.tsv"

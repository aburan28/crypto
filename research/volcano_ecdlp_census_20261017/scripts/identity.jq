def fromhex: ascii_downcase | ltrimstr("0x") | explode | reduce .[] as $c (0; .*16 + (if $c>=97 then $c-87 else $c-48 end));
($tr | map({key:"\(.a),\(.b)", value:.}) | from_entries) as $by
| [ .workloads[] | . as $w | $by["\($w.curve_spec.a),\($w.curve_spec.b)"] as $row
    | ($row.hk[$w.target_index]|split(",")|map(tonumber)) as $hk
    | {ok: (($w.target|map(fromhex)) == $hk), found: ($row != null), gen_ok: ("\($w.curve_spec.gx),\($w.curve_spec.gy)" == $row.g)} ] as $res
| {spec:$n, family:$fam, curves:([.workloads[].curve_spec|"\(.a),\(.b)"]|unique|length), targets:($res|length), mismatches:([$res[]|select((.ok|not) or (.found|not) or (.gen_ok|not))]|length)}

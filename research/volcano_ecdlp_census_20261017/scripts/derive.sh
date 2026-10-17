#!/bin/sh
. "$(dirname "$0")/env.sh"
for V in V1 V2 V3 V4; do
 D=$STUDY/census/$V; mkdir -p $D/transport
 R=$D/walks/D/R.json
 # class table (D walk models)
 jq -s -c --slurpfile r $R --arg V $V '$r[0] as $r | .[] | {volcano:$V, id, level, dist, parent, kernel_x, kernels, nbrs, p:$r.p, a, b, j, order:$r.order, r:$r.n, cofactor:$r.cofactor, trace:$r.trace, ell:$r.ell, endo_order_conductor_ell_valuation:.level, depth:$r.depth, crater_shape:$r.crater, controls:[.controls[]|{u,a,b}]}' $D/walks/D/vertices.jsonl | sed 's/^/ /' | sed 's/^ //' > $D/vertices.jsonl
 cp $D/walks/D/edges.jsonl $D/edges.jsonl
 for W in $D/walks/*/; do
  S=$(basename $W); [ -f $W/R.json ] && [ -f $W/vertices.jsonl ] && [ -f $W/edges.jsonl ] || continue
  jq -s -c --slurpfile r $W/R.json --arg S $S --arg V $V '
   . as $v | $r[0] as $r |
   def path($i): if $v[$i].parent < 0 then [$i] else path($v[$i].parent) + [$i] end;
   .[] | . as $x | path($x.id) as $pa |
   ([range(1; $pa|length)] | map($v[$pa[.]].level - $v[$pa[.]-0].level*0)) as $_ |
   ([range(1; $pa|length)] | map($v[$pa[.]].level - $v[$pa[.-1]].level)) as $dl |
   {volcano:$V, src:$S, id:$x.id, j:$x.j, level:$x.level, a:$x.a, b:$x.b, p:$r.p, order:$r.order, r:$r.n, path_len:$x.dist, path:$pa,
    downs:($dl|map(select(.==1))|length), ups:($dl|map(select(.==-1))|length), horizontals:($dl|map(select(.==0))|length),
    g:$x.g, h:$x.h, hk:($x.hk // null), gn:$x.gn, controls:$x.controls}' $W/vertices.jsonl > $D/transport/$S.jsonl
 done
done

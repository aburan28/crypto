# input: transport rows (slurped) of one family. args: V FAM ts ic maxpass nsrc(array of {level,n}) nD(array)
def curve($c; $a; $b; $gx; $gy): {kind:"prime_explicit",p:$c.p,a:$a,b:$b,group_order:$c.order,r:$c.r,gx:$gx,gy:$gy};
def xy($s): ($s|split(",")|map(tonumber));
. as $rows
| ($rows | map(.level) | unique)[] as $L
| ($rows | map(select(.level==$L))) as $cs
| (([$lvlcount[]|select(.level==$L)|.n]|add // 0) * 8) as $wlround
| ((2000 / $wlround) | ceil) as $rounds
| ( if $FAM=="N" then $cs | map(curve(.; .a; .b; xy(.gn)[0]; xy(.gn)[1]))
    elif $FAM=="ctl" then $cs | map(. as $c | .controls[] | curve($c; .a; .b; xy(.g)[0]; xy(.g)[1]))
    else $cs | map(curve(.; .a; .b; xy(.g)[0]; xy(.g)[1])) end ) as $curves
| (($maxpass / (($rounds+1)*8)) | floor) as $per
| [range(0; ($curves|length); $per)] | to_entries[] | . as $e
| ($curves[$e.value : $e.value+$per]) as $chunk
| ( [ {name:"rho-neg",role:"reference",method:{id:"rho.negation"}},
      {name:"rho-plain",role:"baseline",method:{id:"rho.plain"}},
      {name:"bsgs",role:"candidate",method:{id:"bsgs.textbook"}},
      {name:"bsgs-neg",role:"candidate",method:{id:"bsgs.negation"}},
      {name:"ic",role:"candidate",method:{id:"ic.pipeline",params:{factor_base:$ic.factor_base,oracle:$ic.oracle}}},
      {name:"rho-neg-aa",role:"control",method:{id:"rho.negation"}} ]
    + (if $V=="V3" and $L==0 and ($FAM|IN("D","N","ctl")) then [{name:"rho-cubic",role:"control",method:{id:"rho.cubic"}}] else [] end) ) as $arms
| {name:("\($V)-\($FAM)-L\($L)-c\($e.key)"), V:$V, family:$FAM, level:$L, chunk:$e.key, curves:($chunk|length), targets:8, rounds:$rounds, warmup:1, arms:($arms|length),
   passes:(($chunk|length)*8*($rounds+1)),
   spec:{schema:"ecbench.spec/v1", label:("volcano census \($V) \($FAM) level \($L) chunk \($e.key)"),
         description:"Volcano ECDLP census Phase B (research/volcano_ecdlp_census_20261017/PLAN.md). Known-answer synthetic curves, planted_shared targets; no security claim.",
         workloads:{curves:$chunk,targets_per_curve:8,target_seed:$ts,target_kind:"planted_shared"},
         arms:$arms, measurement:{rounds:$rounds,warmup:1,order:"alternate",seed:20261017,isolation_required:"L2",timeout_seconds:600}}}

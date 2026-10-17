def arms($cfgs): [
 {name:"rho-neg",role:"reference",method:{id:"rho.negation"}},
 {name:"rho-plain",role:"baseline",method:{id:"rho.plain"}},
 {name:"bsgs",role:"candidate",method:{id:"bsgs.textbook"}},
 {name:"bsgs-neg",role:"candidate",method:{id:"bsgs.negation"}}
 ] + ($cfgs | map({name:.id,role:"candidate",method:{id:"ic.pipeline",params:{factor_base:.factor_base,oracle:.oracle}}}))
 + [{name:"rho-neg-aa",role:"control",method:{id:"rho.negation"}}];

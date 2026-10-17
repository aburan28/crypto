import fs from 'node:fs';
const root=process.argv[2],out=root+'/source/order_setup.gp';
if(fs.existsSync(out))throw Error('write-once GP export');
const family=JSON.parse(fs.readFileSync(root+'/data/family.json')),fixtures=JSON.parse(fs.readFileSync(root+'/data/fixtures.json'));
let gp='p='+family.p+';r='+family.subgroup_order+';expectedN='+family.order+';checked=0;\n';
for(const node of family.nodes){
  const cid=node.id;
  gp+='E=ellinit(['+node.a+','+node.b+'],p);tw=getwalltime();tc=getabstime();N=ellcard(E);cpu=getabstime()-tc;wall=getwalltime()-tw;if(N!=expectedN,error("curve order"));print("{\\\"stage\\\":\\\"curve_order\\\",\\\"curve\\\":'+cid+',\\\"order\\\":\\\"",N,"\\\",\\\"cpu_ms\\\":",cpu,",\\\"wall_ms\\\":",wall,"}");\n';
  gp+='tw=getwalltime();tc=getabstime();fac=factor(N);cpu=getabstime()-tc;wall=getwalltime()-tw;if(fac!=[2,2;r,1],error("factorization"));print("{\\\"stage\\\":\\\"factorization\\\",\\\"curve\\\":'+cid+',\\\"cpu_ms\\\":",cpu,",\\\"wall_ms\\\":",wall,"}");\n';
  gp+='tw=getwalltime();tc=getabstime();invariants=ellgroup(E);cpu=getabstime()-tc;wall=getwalltime()-tw;print("{\\\"stage\\\":\\\"group_structure\\\",\\\"curve\\\":'+cid+',\\\"invariants\\\":\\\"",invariants,"\\\",\\\"cpu_ms\\\":",cpu,",\\\"wall_ms\\\":",wall,"}");\n';
  // P is repeated in the workload deliberately, just as in the solver panel.
  for(const row of fixtures[String(cid)])for(const role of ['P','Q']){
    const point='['+row[role].map(x=>'Mod('+x+',p)').join(',')+']';
    gp+='P='+point+';tw=getwalltime();tc=getabstime();o=ellorder(E,P,[N,fac]);cpu=getabstime()-tc;wall=getwalltime()-tw;if(o!=r,error("point order"));checked++;print("{\\\"stage\\\":\\\"point_order\\\",\\\"curve\\\":'+cid+',\\\"instance\\\":'+row.id+',\\\"role\\\":\\\"'+role+'\\\",\\\"returned_order\\\":\\\"",o,"\\\",\\\"cpu_ms\\\":",cpu,",\\\"wall_ms\\\":",wall,"}");\n';
  }
}
gp+='if(checked!=9200,error("coverage"));print("{\\\"stage\\\":\\\"complete\\\",\\\"point_checks\\\":",checked,"}");quit();\n';
// One JavaScript escape layer is removed before GP sees string literals.
gp=gp.replaceAll('\\\\\\"','\\"');fs.writeFileSync(out,gp);console.log(JSON.stringify({script:out,curves:46,point_orders:9200,timer_resolution_ms:1,scalar_witnesses_read:false}));

allocatemem(200000000);
pairs=readstr("research/conductor_gap_cm_orbits_20261015/prepared/pairs.gpvec");
work=readstr("research/conductor_gap_cm_orbits_20261015/inputs/workloads.csv");
ids=List();for(i=2,#work,w=eval(Str("[",work[i],"]"));if(w[13]==0&&!setsearch(Set(Vec(ids)),w[2]),listput(ids,w[2]));if(#ids==10,break));
for(i=1,#ids,id=ids[i];v=eval(pairs[id+1]);l=v[2];p=v[3];N=v[4];T=ellinit([v[5],v[6]],p);for(variant=0,1,K=[0];for(x=1000*(variant+1)+id,1000*(variant+1)+id+1000,if(kronecker(lift(Mod(x,p)^3+v[5]*x+v[6]),p)==1,K=ellmul(T,[Mod(x,p),sqrt(Mod(x,p)^3+v[5]*x+v[6])],5);if(K!=[0]&&ellmul(T,K,l)==[0],break)));if(K==[0]||ellmul(T,K,l)!=[0],error("extra kernel search failed"));key=3000+2*(i-1)+variant;print(key," ",l," ",p," ",v[5]," ",v[6]," ",lift(K[1])," ",lift(K[2])," ",v[9]," ",v[10]," ",v[11]," -1 -1");write("research/conductor_gap_cm_orbits_20261015/prepared/extra-kernel-catalog.jsonl",Str("{\"id\":",key,",\"parent_pair\":",id,",\"D\":-19,\"p\":\"",p,"\",\"N\":\"",N,"\",\"ell\":\"",l,"\",\"a\":\"",v[10],"\",\"b\":\"",v[11],"\",\"Gx\":\"",v[12],"\",\"Gy\":\"",v[13],"\",\"prime_kernel_order\":true}"))));
quit;

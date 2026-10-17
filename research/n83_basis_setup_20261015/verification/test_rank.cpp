#define main verify_main
#include "verify_rank.cpp"
#undef main
// Synthetic algebra-only fixture. Not a performance run or selection input.
std::string arr(const std::vector<uint64_t>&v){std::string s="[";for(auto x:v){if(s.size()>1)s+=",";s+=std::to_string(x);}return s+"]";}
int main(int argc,char**argv){try{
 require(argc==2,"usage: test_rank NEW_DIRECTORY");files::path dir=argv[1];require(!files::exists(dir),"Existing test directory");files::create_directories(dir);
 uint64_t seed=29;std::vector<uint64_t> scalars,seeds;for(int j=0;j<2;j++){seed=seed*6364136223846793005ULL+1442695040888963407ULL;seeds.push_back(seed);scalars.push_back((seed>>11)%(R-1)+1);}
 uint64_t l0=mm(scalars[0],power(5,R-2,R),R), l1=(scalars[1]+R-mm(4,l0,R))%R;std::vector<P> reps{mul(l0,G),mul(l1,G)};
 std::string points="[",labels="[";for(size_t c=0;c<2;c++){P p=reps[c];uint64_t coeff=1;for(int n=0;n<83;n++){for(int sign=0;sign<2;sign++){if(points.size()>1){points+=",";labels+=",";}points+=point_json(sign?neg(p):p);labels+="["+std::to_string(c)+","+std::to_string(sign?R-coeff:coeff)+"]";}p=tau(p);coeff=mm(coeff,LAMBDA,R);}}points+="]";labels+="]";
 std::string base="{\"kind\":\"point_defined_factor_base\",\"n\":83,\"a\":1,\"subgroup_order\":"+std::to_string(R)+",\"field_modulus_low_terms\":[0,2,4,7],\"orbit_columns\":2,\"base_hash\":\"synthetic-audit-only\",\"factor_base_point_coordinates\":"+points+",\"factor_base_point_labels\":"+labels+",\"factor_base_representatives\":["+point_json(reps[0])+","+point_json(reps[1])+"]}\n";
 P target=mul(4,reps[0]), mid=mul(2,reps[0]);std::string witness="\"point_indices\":[0,0,0,0],\"x_codes\":[";for(int j=0;j<4;j++){if(j)witness+=",";witness+=qstr(dec(reps[0].x));}witness+="],\"pinned_intermediates\":["+qstr(dec(mid.x))+","+qstr(dec(mid.x))+"]";
 std::string trace="{\"kind\":\"n83_rank_header\",\"schema\":\"n83-basis-rank/v1\",\"n\":83,\"a\":1,\"r\":"+std::to_string(R)+",\"field_modulus_low_terms\":[0,2,4,7],\"generator\":"+point_json(G)+",\"columns\":2,\"base_hash\":\"synthetic-audit-only\",\"initial_rank_seed\":29}\n";
 for(int j=0;j<2;j++){auto shared="\"attempt\":"+std::to_string(j+1)+",\"scalar\":"+std::to_string(scalars[j])+",\"column\":"+std::to_string(j)+",\"rank_before\":"+std::to_string(j);trace+="{\"kind\":\"n83_rank_attempt_start\","+shared+",\"rank_seed\":"+std::to_string(seeds[j])+"}\n";trace+="{\"kind\":\"n83_rank_attempt_end\","+shared+",\"target\":"+point_json(target)+",\"found\":true,"+witness+",\"row\":"+arr(j?std::vector<uint64_t>{4,1,scalars[j]}:std::vector<uint64_t>{5,0,scalars[j]})+",\"rank_after\":"+std::to_string(j+1)+",\"gained\":true,\"probes\":1,\"states\":1,\"timing_ns\":{\"anchor_construction\":1,\"relation_search_excluding_checks\":1,\"relation_checks\":1,\"elimination\":1}}\n";}
 trace+="{\"kind\":\"n83_rank_solution\",\"rank\":2,\"attempts\":2,\"relations\":2,\"failures\":0,\"probes\":2,\"logs\":"+arr({l0,l1})+"}\n{\"kind\":\"n83_rank_trace_end\"}\n";
 std::string answers="{\"kind\":\"compact_orbit_dlp_target\",\"fixture_index\":0,\"group_verified\":true,\"target\":"+point_json(target)+",\"published_q\":"+point_json(target)+","+witness+",\"recovered_scalar\":"+std::to_string(mm(4,l0,R))+"}\n";
 auto b=dir/"base.jsonl", t=dir/"trace.jsonl", a=dir/"answers.jsonl";create_new(b,base);create_new(t,trace);create_new(a,answers);
 auto run=[&](std::string tag,files::path bb,files::path tt,files::path aa){std::vector<std::string> args{argv[0],bb.string(),tt.string(),aa.string(),"29",(dir/(tag+".receipt.json")).string()};std::vector<char*> ptrs;for(auto&x:args)ptrs.push_back(x.data());return verify_main(6,ptrs.data());};require(run("positive",b,t,a)==0,"Positive fixture rejected");
 size_t tests=0;auto reject=[&](std::string name,std::string input,std::string before,std::string after,char which){auto at=input.find(before);require(at!=std::string::npos,"Mutation anchor missing");input.replace(at,before.size(),after);auto path=dir/(name+".jsonl");create_new(path,input);require(run(name,which=='b'?path:b,which=='t'?path:t,which=='a'?path:a)!=0,"Tamper accepted");tests++;};
 reject("duplicate_key",trace,"\"initial_rank_seed\":29","\"initial_rank_seed\":29,\"initial_rank_seed\":29",'t');
 reject("oversized_seed",trace,"\"initial_rank_seed\":29","\"initial_rank_seed\":340282366920938463463374607431768211485",'t');
 reject("bad_row",trace,"\"row\":[5,0,","\"row\":[6,0,",'t');
 reject("bad_rank",trace,"\"rank_after\":1","\"rank_after\":0",'t');
 reject("bad_seed",trace,"\"initial_rank_seed\":29","\"initial_rank_seed\":30",'t');
 reject("bad_gain",trace,"\"gained\":true","\"gained\":false",'t');
 reject("bad_anchor",trace,"\"column\":0","\"column\":1",'t');
 reject("bad_probe_total",trace,"\"probes\":2","\"probes\":3",'t');
 reject("bad_logs",trace,arr({l0,l1}),arr({(l0+1)%R,l1}),'t');
 reject("bad_label",base,"\"factor_base_point_labels\":[[0,1]","\"factor_base_point_labels\":[[0,2]",'b');
 reject("bad_answer",answers,"\"recovered_scalar\":"+std::to_string(mm(4,l0,R)),"\"recovered_scalar\":1",'a');
 reject("bad_witness",trace,"\"point_indices\":[0,0,0,0]","\"point_indices\":[1,0,0,0]",'t');
 reject("missing_footer",trace,"{\"kind\":\"n83_rank_trace_end\"}\n","",'t');
 reject("truncated",trace,"}\n{\"kind\":\"n83_rank_trace_end\"}\n","}",'t');
 std::cout<<"PASS: one positive and "<<tests<<" negative transcript tests\n";return 0;
 }catch(const std::exception&e){std::cerr<<"TEST FAILURE: "<<e.what()<<"\n";return 1;}}

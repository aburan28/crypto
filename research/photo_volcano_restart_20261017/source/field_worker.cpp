#include "ec.hpp"
#include <filesystem>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/resource.h>
#include <time.h>
extern "C" {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include "../field_backend/include/cryptanalysis/ca_rho.h"
#include "../field_backend/include/cryptanalysis/ca_bsgs.h"
#include "../field_backend/include/cryptanalysis/ca_probe.h"
#pragma GCC diagnostic pop
}
U fcpu(){timespec t{};need(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t)==0,"CPU clock");return U(t.tv_sec)*1000000000+U(t.tv_nsec);}
json fusage(){rusage r{};need(getrusage(RUSAGE_SELF,&r)==0,"rusage");return{{"user_us",r.ru_utime.tv_sec*1000000LL+r.ru_utime.tv_usec},{"system_us",r.ru_stime.tv_sec*1000000LL+r.ru_stime.tv_usec},{"max_rss_kib",r.ru_maxrss},{"minor_faults",r.ru_minflt},{"major_faults",r.ru_majflt},{"voluntary_switches",r.ru_nvcsw},{"involuntary_switches",r.ru_nivcsw}};}
json probe_wire(const ca_probe_counts&c){json j;for(unsigned i=0;i<3;i++){json row;
#define FPROBE(name) row[#name]=c.name[i]
FPROBE(add);FPROBE(sub);FPROBE(mul_plain);FPROBE(mont_mul);FPROBE(mont_square);FPROBE(mont_redc);FPROBE(inverse);FPROBE(euclid_divisions);FPROBE(complement);
#undef FPROBE
j[i==0?"Fp":i==1?"scalar_mod_r":"other_modulus"]=row;}return j;}
void fatomic(const std::string&p,const json&j){need(!std::filesystem::exists(p),"refuse overwrite "+p);save(p+".tmp",j);std::filesystem::rename(p+".tmp",p);}
void fencode(const ca_group&g,ca_elem&out,Point p){U words[4]={p.x,p.y,0,0};need(ca_group_encode(&g,&out,words),"encode");}
json fstats(const ca_stats&s){return{{"group_ops",s.group_ops},{"iterations",s.iterations},{"collisions",s.collisions},{"table_entries",s.table_entries},{"bytes_peak",s.bytes_peak},{"threads",s.threads}};}
int main(int argc,char**argv){try{
need(argc==5,"ROOT rho|bsgs OFFSET STRIDE");std::string root=argv[1],alg=argv[2];int offset=std::stoi(argv[3]),stride=std::stoi(argv[4]);need(offset>=0&&offset<stride&&(alg=="rho"||alg=="bsgs"),"args");auto family=load(root+"/data/family.json"),fixtures=load(root+"/data/fixtures.json"),binding=load(root+"/FIELD_FROZEN.json");int completed=0,failed=0;
for(int cid=0;cid<46;cid++)for(int unit=0;unit<(alg=="rho"?100:10);unit++){
int task=cid*(alg=="rho"?100:10)+unit;if(task%stride!=offset)continue;
if(std::filesystem::exists(root+"/FIELD_STOP")){std::cout<<"FIELD_STOP acknowledged\n";return 0;}
std::string suffix="C"+std::to_string(cid)+(alg=="rho"?"_I":"_G")+std::to_string(unit),job=root+"/field_runs/"+alg+"/"+suffix;std::filesystem::create_directories(job);
int lock=open((job+"/.lock").c_str(),O_CREAT|O_RDWR,0600);need(lock>=0,"job lock");if(flock(lock,LOCK_EX|LOCK_NB)!=0){close(lock);continue;}
if(std::filesystem::exists(job+"/complete.json")){auto c=load(job+"/complete.json");need(c["binding"]==binding&&c["verified"]==true,"field checkpoint binding");close(lock);completed++;continue;}
std::string basepath=root+"/runs/"+alg+"/"+suffix+"/complete.json";need(std::filesystem::exists(basepath),"baseline job incomplete");
std::string attempt=job+"/attempt_"+std::to_string(ns())+"_"+std::to_string(getpid());std::filesystem::create_directories(attempt);fatomic(attempt+"/input.json",{{"curve",cid},{"unit",unit},{"algorithm",alg},{"binding",binding},{"baseline_record_read_before_solver",false}});
auto node=family["nodes"][cid];EC check{val(node["a"]),val(node["b"])};const auto&rows=fixtures[std::to_string(cid)];ca_group g{};need(ca_group_ec_init(&g,PRIME,check.a,check.b,ORDER)==CA_OK,"group");g.cofactor=4;ca_elem base,target;ca_bsgs_table*table=nullptr;json records=json::array();bool all=true;
if(alg=="bsgs"){
fencode(g,base,point(rows[unit*10]["P"]));auto before=fusage();U wall=ns(),cpu=fcpu();ca_probe_begin();int rc=ca_bsgs_table_new(&g,&base,134217728,&table);auto count=ca_probe_end();cpu=fcpu()-cpu;wall=ns()-wall;fatomic(attempt+"/setup.json",{{"status",rc},{"wall_ns",wall},{"process_cpu_ns",cpu},{"counters",probe_wire(count)},{"resources_before",before},{"resources_after",fusage()}});all=rc==CA_OK;
}
const int first=alg=="rho"?unit:unit*10,last=alg=="rho"?unit+1:unit*10+10;
if(all)for(int id=first;id<last;id++){
const auto&row=rows[id];Point P=point(row["P"]),Q=point(row["Q"]);fencode(g,base,P);fencode(g,target,Q);ca_stats stats{};ca_rho_params params{};
if(alg=="rho"){ca_rho_params_default(&params);params.seed=val(row["rho_seed"]);params.threads=1;params.max_ops=2000000000ULL;params.max_table_entries=1000000;}
auto before=fusage();U x=0,wall=ns(),cpu=fcpu();ca_probe_begin();int rc=alg=="rho"?ca_rho_solve(&g,&base,&target,&params,&x,&stats):ca_bsgs_table_solve(table,&target,0,0,&x,&stats);auto count=ca_probe_end();cpu=fcpu()-cpu;wall=ns()-wall;
U verify=ns();bool equation=rc==CA_OK&&check.scalar(x,P)==Q;auto baseline=load(basepath);bool match=false;for(const auto&r:baseline["records"])if(r["instance"]==id){match=r["verified"]==true&&val(r["recovered_k"])==x;const auto&s=r["native_stats"];const auto observed=fstats(stats);for(const auto&[key,value]:observed.items())match=match&&s[key]==value;break;}
json record={{"algorithm",alg},{"curve",cid},{"level",node["level"]},{"instance",id},{"generator",row["generator"]},{"key_index",row["key_index"]},{"status",rc},{"recovered_k",std::to_string(x)},{"native_stats",fstats(stats)},{"instrumented_wall_ns",wall},{"instrumented_process_cpu_ns",cpu},{"counters",probe_wire(count)},{"point_order",std::to_string(ORDER)},{"verified_public_equation",equation},{"matched_baseline_answer_and_all_nonclock_stats",match},{"baseline_read_after_native_return",true},{"verified",equation&&match},{"validation_ns",ns()-verify},{"resources_before",before},{"resources_after",fusage()},{"timing_class","counter-instrumented replay, L0; not baseline timing"}};
fatomic(attempt+"/case_"+std::to_string(id)+".json",record);records.push_back(record);all=all&&equation&&match;std::cout<<"FIELD "<<alg<<" C"<<cid<<" I"<<id<<" matched="<<(equation&&match)<<" group_ops="<<stats.group_ops<<" Fp_M="<<count.mont_mul[0]<<" Fp_S="<<count.mont_square[0]<<" Fp_inv="<<count.inverse[0]<<std::endl;
}
if(table){U wall=ns(),cpu=fcpu();ca_probe_begin();ca_bsgs_table_free(table);auto count=ca_probe_end();cpu=fcpu()-cpu;wall=ns()-wall;fatomic(attempt+"/destruction.json",{{"wall_ns",wall},{"process_cpu_ns",cpu},{"counters",probe_wire(count)}});}
json result={{"binding",binding},{"curve",cid},{"unit",unit},{"algorithm",alg},{"verified",all},{"records",records},{"attempt",std::filesystem::path(attempt).filename().string()}};fatomic(attempt+"/result.json",result);if(all){fatomic(job+"/complete.json",result);completed++;}else failed++;close(lock);
}
std::cout<<"FIELD_FINISHED "<<alg<<" partition="<<offset<<"/"<<stride<<" completed="<<completed<<" failed="<<failed<<std::endl;return failed?2:0;
}catch(const std::exception&e){std::cerr<<"FIELD_ERROR "<<e.what()<<'\n';return 1;}}

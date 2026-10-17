#include "ec.hpp"
extern "C" {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include "../field_backend/include/cryptanalysis/ca_rho.h"
#include "../field_backend/include/cryptanalysis/ca_bsgs.h"
#include "../field_backend/include/cryptanalysis/ca_probe.h"
#pragma GCC diagnostic pop
}
json probe_wire(const ca_probe_counts&c){json j;for(unsigned i=0;i<3;i++){json row;
#define FPROBE(name) row[#name]=c.name[i]
FPROBE(add);FPROBE(sub);FPROBE(mul_plain);FPROBE(mont_mul);FPROBE(mont_square);FPROBE(mont_redc);FPROBE(inverse);FPROBE(euclid_divisions);FPROBE(complement);
#undef FPROBE
j[i==0?"Fp":i==1?"scalar_mod_r":"other_modulus"]=row;}return j;}

int main(int argc,char**argv){try{need(argc==5,"ROOT CURVE INSTANCE rho|bsgs");std::string root=argv[1],alg=argv[4];int cid=std::stoi(argv[2]),id=std::stoi(argv[3]);need(cid>=0&&cid<46&&id>=0&&id<100&&(alg=="rho"||alg=="bsgs"),"arguments");auto family=load(root+"/data/family.json"),fixtures=load(root+"/data/fixtures.json");auto node=family["nodes"][cid],row=fixtures[std::to_string(cid)][id];EC ec{val(node["a"]),val(node["b"])};Point P=point(row["P"]),Q=point(row["Q"]);need(ec.on(P)&&ec.on(Q)&&!P.infinity&&!Q.infinity&&ec.scalar(ORDER,P).infinity&&ec.scalar(ORDER,Q).infinity,"public subgroup inputs");ca_group g{};need(ca_group_ec_init(&g,PRIME,ec.a,ec.b,ORDER)==CA_OK,"curve init");g.cofactor=4;ca_elem base,target;U pw[4]={P.x,P.y,0,0},qw[4]={Q.x,Q.y,0,0};need(ca_group_encode(&g,&base,pw)&&ca_group_encode(&g,&target,qw),"encode");ca_stats st{};U x=0,t=ns();int rc;ca_probe_begin();if(alg=="rho"){ca_rho_params par;ca_rho_params_default(&par);par.seed=val(row["rho_seed"]);par.threads=1;par.max_ops=2000000000ULL;par.max_table_entries=1000000;rc=ca_rho_solve(&g,&base,&target,&par,&x,&st);}else{ca_bsgs_params par;ca_bsgs_params_default(&par);par.table_size=134217728;par.max_ops=2000000000ULL;rc=ca_bsgs_solve(&g,&base,&target,0,0,&par,&x,&st);}auto counted=ca_probe_end();U elapsed=ns()-t;bool verified=rc==0&&ec.scalar(x,P)==Q;json result={{"counter_scope","Source-level modular API calls; counter-overhead timing, not baseline timing"},{"counters",probe_wire(counted)},{"algorithm",alg},{"curve",cid},{"instance",id},{"status",rc},{"verified_public_equation",verified},{"expected_scalar_file_read",false},{"recovered_scalar_exported",false},{"upstream_group_ops",st.group_ops},{"iterations",st.iterations},{"table_entries",st.table_entries},{"table_bytes",st.bytes_peak},{"wall_ns",elapsed},{"bsgs_cost_boundary",alg=="bsgs"?json("independent cold setup+query; not reused-table online timing"):json(nullptr)}};std::cout<<result.dump(2)<<'\n';return verified?0:2;}catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<'\n';return 1;}}

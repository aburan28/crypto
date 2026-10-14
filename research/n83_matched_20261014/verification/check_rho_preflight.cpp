#include "fresh_fixture_helpers.hpp"

inline bool equal_json(const J&a,const J&b){
 if(a.type!=b.type||a.scalar!=b.scalar||a.array.size()!=b.array.size()||a.object.size()!=b.object.size())return false;
 for(size_t i=0;i<a.array.size();++i)if(!equal_json(a.array[i],b.array[i]))return false;
 for(const auto&e:a.object){auto p=b.object.find(e.first);if(p==b.object.end()||!equal_json(e.second,p->second))return false;}return true;
}
inline void verify_pair(const J&old,const J&fresh){
 require(old.at("schema").scalar=="n83-rho-speed-evaluation/v1"&&fresh.at("schema").scalar=="n83-rho-matched-target/v1","Unexpected preflight schemas");
 require(old.at("status").scalar=="PASS_VERIFIED_PUBLIC_TARGET"&&fresh.at("status").scalar=="PASS_VERIFIED_PUBLIC_TARGET"&&truth(old.at("group_verified"))&&truth(fresh.at("group_verified")),"Preflight producer did not verify");
 check_modulus(old);check_modulus(fresh);
 require(num(old.at("n"))==83&&num(old.at("a"))==1&&num(old.at("b"))==1&&num(old.at("lambda"))==2998004823950457ULL,"Old preflight curve/lambda mismatch");
 for(const auto&j:{&old,&fresh}){
  require(j->at("curve_slug").scalar=="icv1-f2m83-t6151469093347-9eba6cd5"&&num(j->at("subgroup_order"))==8569786107849059ULL,"Preflight exact curve/subgroup mismatch");
  require(num(j->at("target_count"))==1&&!truth(j->at("known_answer_input"))&&!truth(j->at("validation_sidecar_opened")),"Preflight is not one public point");
  P g=point(j->at("generator")),q=point(j->at("published_q"));U d=num(j->at("recovered_scalar"));
  require(g.x==decimal("6614244651398184860292015")&&g.y==decimal("8868260410029249391539529")&&!q.inf&&on(q)&&mul(8569786107849059ULL,q).inf,"Preflight point/generator/subgroup invalid");
  require(d>0&&d<8569786107849059ULL&&mul(d,g)==q,"Preflight independent canonical scalar replay failed");
 }
 for(const auto&key:{"published_q","generator","recovered_scalar","walk_seed","jump_seed","optimization_config"})require(equal_json(old.at(key),fresh.at(key)),"Old/new preflight point/scalar/seed/policy mismatch");
 P fixed=point(old.at("published_q"));require(fixed.x==decimal("2276293614448838520573184")&&fixed.y==decimal("278690618123076595657135"),"Preflight point differs from fixed speed-study holdout_1");
 require(old.at("public_input_sha256").scalar=="8cb893707f7eb480d1cd600436005d435f17687027eade3c39aa87ab459ae73d"&&fresh.at("provenance").at("public_input_sha256").scalar==old.at("public_input_sha256").scalar,"Preflight point-only file bytes differ from fixed historical input");
 require(num(old.at("walk_seed"))==202610126241ULL,"Preflight walk seed differs from fixed historical holdout seed");
 const auto&c=fresh.at("optimization_config");require(c.at("backend").scalar=="fold"&&c.at("scalar_kernel").scalar=="reference"&&num(c.at("lanes"))==256&&num(c.at("dp_bits"))==6&&!truth(c.at("escape_fruitless"))&&num(c.at("step_cap_factor"))==2000,"Preflight selected policy mismatch");
 for(const auto&key:{"walk_steps","walks","fruitless_cycles","fruitless_escapes","capped_walks","wasted_merges","table_entries","ideal_steps","automorphism_size"})require(equal_json(old.at("outcome").at(key),fresh.at("outcome").at(key)),"Old/new preflight outcome counter mismatch");
 require(num(old.at("outcome").at("table_payload_lower_bound_bytes"))==49*num(fresh.at("outcome").at("table_entries")),"Preflight table payload lower-bound mismatch");
 require(equal_json(old.at("preparation_charges"),fresh.at("preparation_charges_for_this_seed_only"))&&equal_json(old.at("outcome").at("charges_setup_and_solve"),fresh.at("outcome").at("charges_jump_preparation_and_solve")),"Old/new preflight recorded charge mismatch");
 require(truth(fresh.at("fresh_target_dependent_DP_walk_and_rng_state"))&&!truth(fresh.at("cross_target_collision_sharing")),"Preflight new worker shares target-dependent state");
}
int main(int argc,char**argv){try{
 bool self_test=argc==5&&std::string(argv[1])=="--self-test";require(argc==4||self_test,"usage: check_rho_preflight [--self-test] OLD_RECEIPT NEW_TARGET_JSONL NEW_RECEIPT");
 size_t offset=self_test?2:1;files::path old_file=argv[offset],new_file=argv[offset+1],output=argv[offset+2];require(!files::exists(output),"Refusing existing preflight receipt");
 auto old=read_json(old_file.string());auto raw=lines(new_file);require(raw.size()==1,"Preflight new input must contain exactly one target");auto fresh=raw[0];verify_pair(old,fresh);
 bool counter_rejected=false,scalar_rejected=false;
 if(self_test){
  auto bad=fresh;auto&counter=bad.object.at("outcome").object.at("walk_steps");counter.scalar=dec(num(counter)+1);try{verify_pair(old,bad);}catch(const std::exception&){counter_rejected=true;}
  bad=fresh;bad.object.at("recovered_scalar").scalar="1";try{verify_pair(old,bad);}catch(const std::exception&){scalar_rejected=true;}
  require(counter_rejected&&scalar_rejected,"Preflight checker negative self-test admitted corruption");
 }
 std::ostringstream out;out<<"{\n  \"kind\":"<<qstr(self_test?"native_rho_preflight_negative_self_test":"native_rho_singleton_algorithm_semantic_preflight")<<",\n  \"status\":\"PASS\",\n  \"same_point_scalar_seed_policy\":true,\n  \"common_outcome_fields_checked\":9,\n  \"table_payload_lower_bound_derived_and_checked\":true,\n  \"preparation_and_total_charges_equal\":true,\n  \"independent_canonical_scalar_replays\":2,\n  \"old_receipt_sha256\":"<<qstr(digest(old_file))<<",\n  \"new_target_file_sha256\":"<<qstr(digest(new_file))<<",\n  \"checker_source_sha256\":"<<qstr(digest(__FILE__))<<",\n  \"checker_executable_sha256\":"<<qstr(digest(argv[0]))<<",\n  \"mismatched_counter_rejected\":"<<(self_test?(counter_rejected?"true":"false"):"null")<<",\n  \"mismatched_scalar_rejected\":"<<(self_test?(scalar_rejected?"true":"false"):"null")<<",\n  \"scope\":\"Singleton worker reshape control on one fixed historical public target; equality of exported outcomes and charges, not a full trajectory trace or fresh-panel result\",\n  \"performance_claim\":false\n}\n";create_new(output,out.str());std::cout<<"PASS: same old/new rho public point, scalar, seeds, policy, all9 common outcome fields and preparation/total charges"<<(self_test?"; counter/scalar negative self-tests rejected":"")<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}

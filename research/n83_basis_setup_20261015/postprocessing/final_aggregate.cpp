// Separately versioned aggregation of immutable evidence, including failed arms.
// No solves, rank audits or target generation are repeated here.
#define main sealed_stage_checker_main
#include "../stage_verification/check_stage.cpp"
#undef main
#include "statistics.hpp"

void sealed(const files::path& path,const J& expected) {
 require(expected.type=="string" && digest(path)==expected.scalar,("Evidence hash mismatch: "+path.string()).c_str());
}
bool success(const J&e){return e.at("status").scalar=="COMPLETED_VERIFIED";}
void completion_gate(size_t attempts,size_t complete,size_t answers,size_t ranks,size_t aa){
 require(attempts==55&&complete==50&&answers==92&&ranks==30&&aa==5,
         "Actual cross-stage coverage differs; no publication receipt produced");
}
#ifndef FINAL_AGGREGATE_LIBRARY
int main(int argc,char**argv){try{
 require(argc==3,"usage: final_aggregate PRIVATE_ROOT NEW_OUTPUT");files::path root=argv[1];
 auto select=read_json((root/"training-fixtures/selection_freeze.json").string());
 require(select.at("status").scalar=="TRAINING_REVIEWED"&&select.at("selected_candidate").scalar=="C3","Training selection gate");
 auto cold=read_json((root/"cold-audit/stage_receipt.json").string());
 auto cohort=read_json((root/"cohort-audit/stage_receipt.json").string());
 require(cold.at("status").scalar=="PASS_ALL_DECLARED_STAGE_EVIDENCE"&&num(cold.at("completed"))==26,"Cold gate");
 require(cohort.at("status").scalar=="FAILED_EVIDENCE_RETAINED"&&num(cohort.at("completed"))==8&&num(cohort.at("invalid_or_missing"))==3,"Retained incomplete cohort gate");
 std::map<std::string,std::string> cost_seals;
 for(auto*receipt:{&select,&cold,&cohort})for(auto&v:receipt->at("process_costs").array)cost_seals[v.at("id").scalar]=v.at("execution_sha256").scalar;
 std::ostringstream runs,profiles,aa_rows,bindings;runs<<"[";profiles<<"[";aa_rows<<"[";bindings<<"[";
 bool fr=true,fp=true,fa=true,fb=true;size_t attempts=0,complete=0,answers=0,ic_answers=0,rho_answers=0,ranks=0,aa_count=0;
 std::map<std::string,J> aa_work;std::map<std::string,double> aa_time;
 std::map<std::string,std::vector<double>> wall,cpus,rss,online;
 std::array<std::vector<double>,3> cohort_costs;
 std::array<std::array<std::vector<double>,8>,2> target_times;
 auto bind=[&](const files::path&p){if(!fb)bindings<<",";fb=false;bindings<<"{\"path\":"<<qstr(p.string())<<",\"sha256\":"<<qstr(digest(p))<<"}";};
 for(auto p:{root/"training-fixtures/selection_freeze.json",root/"cold-audit/stage_receipt.json",root/"cohort-audit/stage_receipt.json",root/"cold_statistics.json",root/"training_profiles.json"})bind(p);
 for(auto stage:{std::string("training"),std::string("cold"),std::string("cohort")}){
  auto campaign=root/(stage+"-campaign");auto matrix_path=stage=="training"?root/"training-campaign/matrix.json":root/"holdout-plan/matrix.json";
  auto matrix=read_json(matrix_path.string());bind(matrix_path);require(digest(campaign/"matrix.json")==digest(matrix_path),"Campaign matrix mismatch");
  const auto&receipt=stage=="training"?select:stage=="cold"?cold:cohort;sealed(matrix_path,receipt.at("matrix_sha256"));
  for(auto&[p,h]:matrix.at("input_bindings").object){sealed(p,h);bind(p);}
  sealed(matrix.at("config").at("ic_binary").scalar,matrix.at("ic_binary_sha256"));
  sealed(matrix.at("config").at("rho_binary").scalar,matrix.at("rho_binary_sha256"));
  for(auto&row:matrix.at("rows").array){if(stage!="training"&&row.at("stage").scalar!=stage)continue;
   auto id=row.at("id").scalar;auto dir=campaign/id;auto e=read_json((dir/"execution.json").string());bool ic=row.at("method").scalar=="ic",ok=success(e),is_aa=id.rfind("aa_",0)==0;
   envelope(e);bind(dir/"execution.json");sealed(dir/"stdout.jsonl",e.at("stdout_sha256"));sealed(dir/"stderr.txt",e.at("stderr_sha256"));
   require(e.at("binary_sha256").scalar==matrix.at(ic?"ic_binary_sha256":"rho_binary_sha256").scalar,"Worker identity mismatch");
   if(cost_seals.count(id))require(cost_seals[id]==digest(dir/"execution.json"),"Audit cost seal mismatch");
   else require(!ok&&stage=="cohort"&&!ic,"Unbound completed cost");
   std::string role=stage=="training"?(ic?row.at("candidate").at("id").scalar:"rho"):row.at("role").scalar;
   std::string group=stage+"/"+role;attempts++;if(ok)complete++;
   if(!fr)runs<<",";fr=false;runs<<"{\"stage\":"<<qstr(stage)<<",\"id\":"<<qstr(id)<<",\"role\":"<<qstr(role)<<",\"completed\":"<<(ok?"true":"false")<<",\"wall_s\":"<<e.at("elapsed_process_wall_s").scalar<<",\"cpu_s\":"<<e.at("cpu_total_s").scalar<<",\"peak_rss_bytes\":"<<e.at("peak_rss_bytes").scalar<<",\"stderr\":"<<qstr(bytes(dir/"stderr.txt"))<<",\"execution_sha256\":"<<qstr(digest(dir/"execution.json"))<<"}";
   if(!ok){require((stage=="training"&&role=="C4")||(stage=="cohort"&&!ic),"Unexpected failure");continue;}
   sealed(dir/"answers.jsonl",e.at("semantic_sha256"));auto answer_rows=lines(dir/"answers.jsonl");require(answer_rows.size()==num(row.at("target_count")),"Answer count");answers+=answer_rows.size();if(ic)ic_answers+=answer_rows.size();else rho_answers+=answer_rows.size();
   if(ic){auto audit_path=stage=="training"?root/"training-audit"/(id+".json"):root/(stage+"-audit")/(id+".rank.json");auto audit=read_json(audit_path.string());require(audit.at("status").scalar=="PASS"&&num(audit.at("answers"))==answer_rows.size(),"Rank audit gate");sealed(dir/"base.jsonl",audit.at("base_sha256"));sealed(dir/"rank.jsonl",audit.at("trace_sha256"));sealed(dir/"answers.jsonl",audit.at("answers_sha256"));bind(audit_path);ranks++;}
   if(is_aa){auto pair=id.substr(0,id.size()-2);auto w=real(e.at("elapsed_process_wall_s"));if(!aa_work.count(pair)){aa_work[pair]=answer_rows[0].at("outcome");aa_time[pair]=w;}else{require(same(aa_work[pair],answer_rows[0].at("outcome")),"Recursive A/A work mismatch");if(!fa)aa_rows<<",";fa=false;aa_rows<<std::setprecision(17)<<"{\"pair\":"<<qstr(pair)<<",\"stage\":"<<qstr(stage)<<",\"first_wall_s\":"<<aa_time[pair]<<",\"second_wall_s\":"<<w<<",\"larger_over_smaller\":"<<std::max(w,aa_time[pair])/std::min(w,aa_time[pair])<<",\"identical_recursive_work\":true}";aa_count++;}continue;}
   wall[group].push_back(real(e.at("elapsed_process_wall_s")));cpus[group].push_back(real(e.at("cpu_total_s")));rss[group].push_back(real(e.at("peak_rss_bytes")));
   if(stage=="cohort"){require(ic,"No cohort rho answers expected");cohort_costs[role=="control"?0:1].push_back(real(e.at("elapsed_process_wall_s")));}
   size_t clock_index=0;auto stdout_rows=lines(dir/"stdout.jsonl");for(auto&v:stdout_rows)if(v.at("kind").scalar=="n83_aligned_query_receipt"){double t=double(num(v.at("interval_ns")))/1e9;online[group].push_back(t);if(stage=="cohort"){auto target=size_t(num(row.at("target_indices").at(clock_index)));target_times[role=="control"?0:1][target].push_back(t);}clock_index++;}
   require(clock_index==answer_rows.size(),"Aligned clock coverage mismatch");
   if(!ic)continue;
   if(!fp)profiles<<",";fp=false;profiles<<"{\"stage\":"<<qstr(stage)<<",\"id\":"<<qstr(id)<<",\"role\":"<<qstr(role)<<",\"summary\":";
   bool found=false;for(auto&v:stdout_rows)if(v.at("kind").scalar=="compact_orbit_dlp_summary"){profiles<<dump(v);found=true;}require(found,"Missing producer summary");
   std::map<std::string,U> phases;U probes=0,states=0,gain=0,useful=0,attempt_count=0;profiles<<",\"attempt_rows\":[";bool first=true;
   for(auto&v:lines(dir/"rank.jsonl")){auto k=v.at("kind").scalar;if(k=="n83_rank_attempt_end"){attempt_count++;probes+=num(v.at("probes"));states+=num(v.at("states"));gain+=truth(v.at("gained"));useful+=truth(v.at("found"));for(auto&[key,t]:v.at("timing_ns").object)phases[key]+=num(t);if(!first)profiles<<",";first=false;profiles<<"{";bool fk=true;for(auto key:{"attempt","found","gained","probes","states","rank_before","rank_after","timing_ns","memory"}){if(!fk)profiles<<",";fk=false;profiles<<qstr(key)<<":"<<dump(v.at(key));}profiles<<"}";}else if(k=="n83_rank_solution")phases["backsolve"]+=num(v.at("backsolve_ns"));else if(k=="n83_rank_trace_end")phases["trace_prior_format_write_flush"]+=num(v.at("prior_trace_write_flush_ns"));}
   profiles<<"],\"attempts\":"<<dec(attempt_count)<<",\"probes\":"<<dec(probes)<<",\"states\":"<<dec(states)<<",\"rank_gain\":"<<dec(gain)<<",\"useful_relations\":"<<dec(useful)<<",\"phase_sum_ns\":{";bool first_phase=true;for(auto&[k,v]:phases){if(!first_phase)profiles<<",";first_phase=false;profiles<<qstr(k)<<":"<<dec(v);}profiles<<"}}";
  }
 }
 completion_gate(attempts,complete,answers,ranks,aa_count);require(ic_answers==72&&rho_answers==20,"Answer method split mismatch");
 cohort_costs[2]=cohort_costs[0];auto comparison=paired(cohort_costs,202610157802ULL)[2];
 std::ostringstream out;out<<std::setprecision(17)<<"{\"kind\":\"n83_basis_final_evidence/v2-incomplete-cohort\",\"status\":\"COMPLETE_ACCOUNTING_WITH_UNRESOLVED_COHORT_REFERENCE\",\"planned_processes\":55,\"planned_answers\":118,\"attempted\":"<<attempts<<",\"completed\":"<<complete<<",\"failed\":5,\"canonical_answers\":"<<answers<<",\"ic_answers\":"<<ic_answers<<",\"rho_answers\":"<<rho_answers<<",\"rank_log_audits\":"<<ranks<<",\"AA_pairs\":"<<aa_count<<",\"cohort_IC_over_rho\":null,\"cohort_selected_over_control\":{\"mean_cost_ratio\":"<<comparison.ratio<<",\"conditional_95_percent_interval\":["<<comparison.lo<<","<<comparison.hi<<"],\"bootstrap_draws\":20000,\"bootstrap_seed\":202610157802},\"summaries\":{";
 bool first_summary=true;for(auto&[k,v]:wall){if(!first_summary)out<<",";first_summary=false;out<<qstr(k)<<":{\"wall_s\":"<<descriptive(v)<<",\"cpu_s\":"<<descriptive(cpus[k])<<",\"peak_rss_bytes\":"<<descriptive(rss[k])<<",\"aligned_online_s\":"<<descriptive(online[k])<<"}";}
 out<<"},\"online_cold_selected_over_rho\":"<<mean(online["cold/selected"])/mean(online["cold/rho"])<<",\"cohort_online_target_rows\":[";
 for(size_t target=0;target<8;target++){if(target)out<<",";out<<"{\"target_index\":"<<target<<",\"rho\":null";for(size_t role=0;role<2;role++){auto&v=target_times[role][target];require(v.size()==3,"Three repeat IC online coverage");out<<","<<qstr(role==0?"control":"selected")<<":{\"mean_s\":"<<mean(v)<<",\"recorded_repeats_s\":[";for(size_t i=0;i<v.size();i++){if(i)out<<",";out<<v[i];}out<<"]}";}out<<"}";}
 out<<"],\"all_process_rows\":"<<runs.str()<<"],\"profiles\":"<<profiles.str()<<"],\"AA_spreads\":"<<aa_rows.str()<<"],\"bindings\":"<<bindings.str()<<"],\"source_sha256\":"<<qstr(digest(__FILE__))<<",\"binary_sha256\":"<<qstr(digest(argv[0]))<<"}\n";
 create_new(argv[2],out.str());std::cout<<"Bound 55 attempts / 50 complete / 92 answers / 30 rank audits / 5 A/A pairs; cohort rho unavailable\n";return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
#endif

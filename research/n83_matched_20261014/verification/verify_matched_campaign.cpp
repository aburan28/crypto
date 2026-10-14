// Independent arithmetic and retained-file verification. No solver is invoked.
#include "matched_fixture_helpers.hpp"

inline bool equal_json(const J&a,const J&b) {
 if(a.type!=b.type||a.scalar!=b.scalar||a.array.size()!=b.array.size()||a.object.size()!=b.object.size())return false;
 for(size_t i=0;i<a.array.size();++i)if(!equal_json(a.array[i],b.array[i]))return false;
 for(const auto&e:a.object){auto p=b.object.find(e.first);if(p==b.object.end()||!equal_json(e.second,p->second))return false;}return true;
}
struct Run{std::string id,kind,arm,leg;size_t cohort=0,pair=0;std::vector<std::string> targets;std::vector<uint64_t> seeds;};
inline std::string run_id(const std::string&prefix,size_t n,const std::string&suffix){std::ostringstream s;s<<prefix<<std::setfill('0')<<std::setw(2)<<n<<"-"<<suffix;return s.str();}
inline std::vector<Run> independent_matrix(){
 std::vector<Run> out;
 for(size_t pair=1;pair<=5;++pair)for(const auto&leg:{"A","B"})out.push_back({run_id("AA",pair,leg),"aa","rho",leg,0,pair,{matched_id(0)},{202610146410ULL+pair}});
 for(size_t cohort=1;cohort<=3;++cohort){
  std::vector<size_t> indices;for(size_t i=0;i<32;++i)indices.push_back((i+(cohort-1)*11)%32);if(cohort==2)std::reverse(indices.begin(),indices.end());
  std::vector<std::string> ids;std::vector<uint64_t> seeds;for(auto i:indices){ids.push_back(matched_id(i));seeds.push_back(202610146501ULL+(cohort-1)*32+i);}
  std::vector<std::string> arms=cohort%2?std::vector<std::string>{"rho","ic"}:std::vector<std::string>{"ic","rho"};
  for(const auto&arm:arms)out.push_back({run_id("C",cohort,arm),"cohort",arm,"",cohort,0,ids,seeds});
  size_t first=cohort==1?1:(cohort==2?4:7),last=cohort==1?3:(cohort==2?6:8);
  for(size_t pair=first;pair<=last;++pair){arms=pair%2?std::vector<std::string>{"ic","rho"}:std::vector<std::string>{"rho","ic"};for(const auto&arm:arms)out.push_back({run_id("K",pair,arm),"cold",arm,"",0,pair,{matched_id(pair-1)},{202610146600ULL+pair}});}
 }require(out.size()==32,"Independent run matrix size failure");return out;
}
inline double finite_nonnegative(const J&j,const std::string&key){double x=real(j.at(key));require(std::isfinite(x)&&x>=0,"Nonfinite/negative time field");return x;}
struct CheckedRow{
 std::string run,id,arm,status,error,scalar;bool present=false,answer=false,replay=false,witness=false,sidecar=false,complete=false;
 double online=0,phase_sum=0,gap=0;
};

int main(int argc,char**argv){try{
 require(argc==8,"usage: verify_matched_campaign ARCHIVE PRIOR_OPTIMIZATION_FIXTURES PRIOR_SPEED_FIXTURES MATCHED_FIXTURES CAMPAIGN FROZEN_PROTOCOL NEW_RECEIPT");
 files::path archive=argv[1],prior=argv[2],speed=argv[3],fresh=argv[4],campaign=argv[5],protocol_file=argv[6],output=argv[7];
 require(!files::exists(output),"Refusing existing final verification receipt");
 auto instance=verify_instance(archive);verify_and_exclude_both_prior_panels(prior,speed,instance);auto fixtures=verify_matched_fixtures(fresh,prior,speed,instance);
 auto protocol=read_json(protocol_file.string()),manifest=read_json((campaign/"campaign_manifest.json").string());auto manifest_file=campaign/"campaign_manifest.json";
 require(protocol.at("schema").scalar=="n83.matched-campaign/v1"&&manifest.at("schema").scalar=="n83.matched-manifest/v1"&&manifest.at("protocol_sha256").scalar==digest(protocol_file)&&equal_json(manifest.at("protocol"),protocol),"Frozen protocol/manifest binding failure");
 require(num(manifest.at("requested_child_processes"))==32&&num(manifest.at("requested_target_answers"))==218&&num(manifest.at("target_count"))==32&&num(manifest.at("shared_setup_repetitions"))==3&&num(manifest.at("true_cold_target_subset_count"))==8&&num(manifest.at("aa_pairs"))==5,"Frozen requested matrix counts mismatch");
 const auto matrix=independent_matrix();const auto&ordered=manifest.at("ordered_runs").array;require(ordered.size()==32,"Frozen ordered run count mismatch");
 files::path cwd=protocol.at("cwd").scalar,rho_binary=protocol.at("rho_binary").scalar,ic_binary=protocol.at("ic_binary").scalar,base_file=protocol.at("base_file").scalar;
 require(digest(base_file)==digest(archive/"base_n83_K600.jsonl"),"IC base differs from independently validated K600");
 std::map<files::path,std::string> snapshots;std::set<files::path> bound;size_t hash_checks=0,source_checks=0,env_checks=0,resource_checks=0,observation_checks=0;
 auto hash=[&](const files::path&p){auto canonical=files::weakly_canonical(p);auto f=snapshots.find(canonical);if(f!=snapshots.end())return f->second;auto h=digest(p);snapshots[canonical]=h;return h;};
 auto match=[&](const files::path&p,const std::string&expected){require(hash(p)==expected,"Recorded file SHA256 binding mismatch");bound.insert(files::weakly_canonical(p));++hash_checks;};
 match(protocol_file,manifest.at("protocol_sha256").scalar);
 for(const auto&pair:std::vector<std::pair<std::string,std::string>>{{"rho_binary","expected_rho_binary_sha256"},{"ic_binary","expected_ic_binary_sha256"},{"base_file","expected_base_file_sha256"},{"validation_sidecar","expected_validation_sidecar_sha256"},{"generator_receipt","expected_generator_receipt_sha256"},{"human_protocol","human_protocol_sha256"}})match(protocol.at(pair.first).scalar,protocol.at(pair.second).scalar);
 for(const auto&e:manifest.at("generated_input_sha256").object)match(e.first,e.second.scalar);
 const auto&provenance=manifest.at("provenance");for(const auto&e:provenance.at("input_binary_bindings").object)match(e.first,e.second.scalar);match(provenance.at("controller_executable").scalar,provenance.at("controller_executable_sha256").scalar);
 auto runner=cwd/"research/n83_matched_20261014/native_runner";
 for(const auto&pair:std::vector<std::pair<std::string,std::string>>{{"main","src/main.rs"},{"monitor","src/monitor.rs"},{"manifest","Cargo.toml"},{"lock","Cargo.lock"}}){match(runner/pair.second,provenance.at("compiled_source_sha256").at(pair.first).scalar);++source_checks;}
 if(has_field(protocol,"relevant_source_sha256"))for(const auto&e:protocol.at("relevant_source_sha256").object){files::path p=e.first;match(p.is_absolute()?p:cwd/p,e.second.scalar);++source_checks;}
 if(has_field(protocol,"build_receipt")&&has_field(protocol,"build_receipt_sha256"))match(protocol.at("build_receipt").scalar,protocol.at("build_receipt_sha256").scalar);
 auto config_file=campaign/"inputs/rho_config.json";auto config=read_json(config_file.string());
 require(config.at("backend").scalar=="fold"&&config.at("scalar_kernel").scalar=="reference"&&num(config.at("lanes"))==256&&num(config.at("dp_bits"))==6&&!truth(config.at("escape_fruitless"))&&num(config.at("step_cap_factor"))==2000&&equal_json(config,manifest.at("rho_policy")),"Frozen selected rho policy mismatch");
 for(const auto&f:fixtures)match(fresh/(f.id+".jsonl"),protocol.at("expected_fixture_file_sha256").at(f.id).scalar);
 const std::map<std::string,std::string> clean_env={{"KIC_RANK_THREADS","1"},{"KIC_RHO_THREADS","1"},{"KIC_TARGET_THREADS","1"},{"KIC_THREADS","1"},{"LC_ALL","C"},{"MKL_NUM_THREADS","1"},{"OMP_NUM_THREADS","1"},{"OPENBLAS_NUM_THREADS","1"},{"PATH","/usr/bin:/bin"},{"RAYON_NUM_THREADS","1"}};
 std::vector<CheckedRow> rows;std::vector<std::string> errors;std::set<files::path> planned_dirs;std::map<std::string,J> aa_records;
 size_t processes=0,completed=0,checked=0,missing=0,invalid=0,nonanswers=0,ic_answers=0,rho_answers=0,witnesses=0;
 U previous_finished=0;
 for(size_t ri=0;ri<matrix.size();++ri){
  const auto&r=matrix[ri];const auto&declared=ordered[ri];auto dir=campaign/"runs"/r.id,point_file=campaign/"inputs"/(r.id+"-points.jsonl"),seed_file=campaign/"inputs"/(r.id+"-seeds.json");planned_dirs.insert(files::weakly_canonical(dir));
  require(declared.at("run_id").scalar==r.id&&declared.at("kind").scalar==r.kind&&declared.at("arm").scalar==r.arm&&declared.at("leg").scalar==r.leg&&num(declared.at("cohort_id"))==r.cohort&&num(declared.at("pair_id"))==r.pair&&files::path(declared.at("out_dir").scalar)==dir&&files::path(declared.at("target_file").scalar)==point_file,"Frozen run order/identity/path mismatch");
  const auto&ids=declared.at("target_ids").array;auto input=lines(point_file);require(ids.size()==r.targets.size()&&input.size()==r.targets.size(),"Frozen target/input list count mismatch");
  for(size_t i=0;i<r.targets.size();++i){require(ids[i].scalar==r.targets[i],"Frozen target order mismatch");auto f=std::find_if(fixtures.begin(),fixtures.end(),[&](const Fixture&f){return f.id==r.targets[i];});require(f!=fixtures.end()&&point(input[i])==f->point&&input[i].type=="array"&&input[i].array.size()==2,"Frozen matrix point differs from validated fixture");}
  if(r.arm=="rho"){require(files::path(declared.at("seed_file").scalar)==seed_file&&declared.at("rank_seed").scalar=="null","Frozen rho seed path/rank input mismatch");auto s=read_json(seed_file.string());require(s.array.size()==r.seeds.size()&&declared.at("walk_seeds").array.size()==r.seeds.size(),"Frozen rho seed count mismatch");for(size_t i=0;i<r.seeds.size();++i)require(num(s.at(i))==r.seeds[i]&&num(declared.at("walk_seeds").at(i))==r.seeds[i],"Frozen rho seed order/value mismatch");}
  else require(num(declared.at("rank_seed"))==29&&declared.at("walk_seeds").scalar=="null"&&declared.at("seed_file").scalar=="null","Frozen IC rank seed mismatch");
  auto execution_file=dir/"execution.json";
  if(!files::is_regular_file(execution_file)){errors.push_back("Missing execution: "+r.id);for(const auto&id:r.targets){CheckedRow row;row.run=r.id;row.id=id;row.arm=r.arm;row.status="NOT_ATTEMPTED_OR_MISSING_RECEIPT";rows.push_back(row);++missing;}continue;}
  ++processes;auto execution=read_json(execution_file.string());auto status=execution.at("status").scalar;
  bool process_ok=status=="COMPLETED_VERIFIED"&&execution.at("return_code").scalar=="0"&&execution.at("termination_signal").scalar=="null"&&execution.at("timeout").scalar=="false";if(process_ok)++completed;
  auto targets_file=dir/"targets.jsonl";auto semantic=files::is_regular_file(targets_file)?lines(targets_file):std::vector<J>{};std::string metadata_error;J observation;
  try{
   require(execution.at("schema").scalar=="n83.process-receipt/v1"&&!truth(execution.at("known_answer_environment_present")),"Process schema/known-answer environment mismatch");
   const auto&args=execution.at("argv").array;require(args.size()==5,"Unexpected solver argument count");auto binary=r.arm=="rho"?rho_binary:ic_binary;
   require(files::path(args[0].scalar)==binary&&files::path(args[4].scalar)==targets_file,"Solver executable/output differs from matrix");
   for(const auto&a:args)require(a.scalar.find("validation_sidecar")==std::string::npos,"Validation scalar sidecar supplied to measured child");
   if(r.arm=="rho")require(files::path(args[1].scalar)==point_file&&files::path(args[2].scalar)==seed_file&&files::path(args[3].scalar)==config_file,"Rho argv differs from matrix");else require(files::path(args[1].scalar)==base_file&&files::path(args[2].scalar)==point_file&&args[3].scalar=="29","IC argv differs from matrix");
   match(binary,execution.at("binary_sha256").scalar);require(execution.at("binary_sha256").scalar==protocol.at(r.arm=="rho"?"expected_rho_binary_sha256":"expected_ic_binary_sha256").scalar,"Executed binary differs from protocol");
   match(point_file,execution.at("target_sha256").scalar);for(const auto&e:execution.at("argv_existing_file_sha256").object)match(e.first,e.second.scalar);
   require(execution.at("argv_existing_file_sha256").object.count(point_file.string())&&execution.at("argv_existing_file_sha256").object.count(r.arm=="rho"?config_file.string():base_file.string()),"Missing point/config/base launch hash");
   require(files::path(execution.at("semantic_path").scalar)==targets_file&&files::path(execution.at("stdout_path").scalar)==dir/"stdout.jsonl"&&files::path(execution.at("stderr_path").scalar)==dir/"stderr.txt","Process output path mismatch");
   match(dir/"stdout.jsonl",execution.at("stdout_sha256").scalar);match(dir/"stderr.txt",execution.at("stderr_sha256").scalar);if(files::is_regular_file(targets_file))match(targets_file,execution.at("semantic_sha256").scalar);
   std::map<std::string,std::string> env;for(const auto&e:execution.at("environment").object)env[e.first]=e.second.scalar;require(env==clean_env,"Recorded clean child environment mismatch");++env_checks;
   const auto&envp=execution.at("resource_envelope");require(num(envp.at("cpu_count"))==1&&envp.at("cpu_affinity").array.size()==1&&num(envp.at("cpu_affinity").at(0))==3&&num(envp.at("threads"))==1&&num(envp.at("rank_threads"))==1&&num(envp.at("online_target_threads"))==1&&num(envp.at("address_space_limit_bytes"))==12884901888ULL&&real(envp.at("timeout_s"))==2100.,"CPU/thread/memory/timeout envelope mismatch");++resource_checks;
   finite_nonnegative(execution,"elapsed_process_wall_s");require(std::abs(finite_nonnegative(execution,"cpu_user_s")+finite_nonnegative(execution,"cpu_system_s")-finite_nonnegative(execution,"cpu_total_s"))<1e-9&&num(execution.at("peak_rss_native_kib"))*1024==num(execution.at("peak_rss_bytes")),"CPU sum/RSS conversion mismatch");
   require(num(execution.at("started_epoch_ms"))>=num(manifest.at("frozen_before_any_child_epoch_ms"))&&num(execution.at("started_epoch_ms"))>=previous_finished&&num(execution.at("finished_epoch_ms"))>=num(execution.at("started_epoch_ms")),"Actual run order overlaps or predates manifest");previous_finished=num(execution.at("finished_epoch_ms"));
   auto started=read_json((dir/"started.json").string());for(const auto&key:{"argv","cwd","environment","resource_envelope","target_sha256","binary_sha256","argv_existing_file_sha256","started_epoch_ms"})require(equal_json(started.at(key),execution.at(key)),"Started/execution metadata mismatch");
   observation=read_json((dir/"observation.json").string());match(execution_file,observation.at("execution_receipt_sha256").scalar);require(observation.at("run_id").scalar==r.id&&observation.at("kind").scalar==r.kind&&observation.at("arm").scalar==r.arm&&observation.at("campaign_manifest_sha256").scalar==hash(manifest_file)&&equal_json(observation.at("raw_process_receipt"),execution)&&observation.at("target_rows").array.size()==r.targets.size(),"Observation/manifest/execution mismatch");
   for(const auto&key:{"elapsed_process_wall_s","cpu_user_s","cpu_system_s","cpu_total_s","peak_rss_bytes"})require(equal_json(observation.at(key),execution.at(key)),"Observation process cost differs from execution receipt");
   require(num(observation.at("target_count"))==r.targets.size()&&num(observation.at("cohort_id"))==r.cohort&&num(observation.at("pair_id"))==r.pair&&observation.at("leg").scalar==r.leg&&equal_json(observation.at("target_ids"),declared.at("target_ids")),"Observation matrix metadata mismatch");++observation_checks;
   auto stdout_rows=lines(dir/"stdout.jsonl");std::vector<J> summaries;for(const auto&s:stdout_rows)if(has_field(s,"kind")&&s.at("kind").scalar==(r.arm=="rho"?"n83_rho_matched_summary":"compact_orbit_dlp_summary"))summaries.push_back(s);
   if(process_ok){
    require(summaries.size()==1&&semantic.size()==r.targets.size(),"Completed process lacks exact summary/target count");const auto&s=summaries[0];require(num(s.at("targets"))==r.targets.size()&&num(s.at("targets_solved"))==r.targets.size()&&num(s.at("targets_failed"))==0&&equal_json(observation.at("summary"),s),"Completed summary counts mismatch");
    if(r.arm=="ic"){require(num(s.at("n"))==83&&num(s.at("a"))==1&&num(s.at("rank"))==600&&num(s.at("orbit_columns"))==600&&num(s.at("factor_base_points"))==99600&&s.at("base_hash").scalar==protocol.at("expected_base_hash").scalar&&num(s.at("threads"))==1&&num(s.at("target_threads"))==1,"IC summary curve/base/rank/thread mismatch");for(const auto&key:{"setup_total","rank_stage","linear_algebra","targets_total","process_total"})finite_nonnegative(s.at("timing_ms"),key);}
    else{require(equal_json(s.at("optimization_config"),config)&&!truth(s.at("target_dependent_state_shared"))&&num(s.at("target_independent_curve_normal_basis_constructions"))==1&&num(s.at("G_only_jump_tables"))==r.targets.size(),"Rho summary policy/shared state mismatch");finite_nonnegative(s,"producer_preparation_ms");}
    double preparation=finite_nonnegative(s.at("timing_ms"),"setup_total");if(r.arm=="ic"){preparation+=finite_nonnegative(s.at("timing_ms"),"rank_stage");preparation+=finite_nonnegative(s.at("timing_ms"),"linear_algebra");}
    require(std::abs(finite_nonnegative(observation,"producer_preparation_ms")-preparation)<1e-6&&observation.at("status").scalar=="COMPLETED_VERIFIED"&&truth(observation.at("all_verified")),"Observation preparation cost/completion status mismatch");
   }
  }catch(const std::exception&e){metadata_error=e.what();errors.push_back(r.id+": "+metadata_error);++invalid;}
  if(semantic.size()!=r.targets.size())errors.push_back(r.id+": retained target count "+std::to_string(semantic.size())+" differs from requested "+std::to_string(r.targets.size()));
  for(size_t i=0;i<r.targets.size();++i){
   auto f=std::find_if(fixtures.begin(),fixtures.end(),[&](const Fixture&f){return f.id==r.targets[i];});CheckedRow row;row.run=r.id;row.id=r.targets[i];row.arm=r.arm;row.status=status;
   if(i>=semantic.size()){++missing;rows.push_back(row);continue;}row.present=true;++checked;const auto&raw=semantic[i];
   try{
    require(point(raw.at("published_q"))==f->point&&point(raw.at("target"))==f->point&&point(raw.at("generator"))==instance.g&&num(raw.at("fixture_index"))==i,"Target row order/point/generator/index mismatch");
    if(!metadata_error.empty())throw std::runtime_error("Recorded process metadata invalid: "+metadata_error);
    const auto&observed=observation.at("target_rows").at(i);
    require(observed.at("run_id").scalar==r.id&&observed.at("target_id").scalar==r.targets[i]&&observed.at("arm").scalar==r.arm&&observed.at("kind").scalar==r.kind&&num(observed.at("cohort_id"))==r.cohort&&num(observed.at("pair_id"))==r.pair&&observed.at("leg").scalar==r.leg&&equal_json(observed.at("rawtargetrecord"),raw),"Observation target identity or raw record differs from producer");
    if(r.arm=="ic")require(num(raw.at("n"))==83&&num(raw.at("a"))==1&&raw.at("published_fixture_scalar").scalar=="null","IC target identity/known scalar mismatch");
    else{
     check_modulus(raw);require(raw.at("schema").scalar=="n83-rho-matched-target/v1"&&num(raw.at("subgroup_order"))==instance.r&&num(raw.at("target_count"))==1&&num(raw.at("cohort_target_count"))==r.targets.size()&&num(raw.at("walk_seed"))==r.seeds[i]&&equal_json(raw.at("optimization_config"),config)&&!truth(raw.at("known_answer_input"))&&!truth(raw.at("validation_sidecar_opened"))&&truth(raw.at("fresh_target_dependent_DP_walk_and_rng_state"))&&!truth(raw.at("cross_target_collision_sharing")),"Rho target seed/policy/input/state reset mismatch");
     const auto&p=raw.at("provenance");match(point_file,p.at("public_input_sha256").scalar);match(seed_file,p.at("seed_input_sha256").scalar);match(config_file,p.at("config_input_sha256").scalar);match(rho_binary,p.at("executable_sha256").scalar);require(p.at("compiled_source_sha256").object.size()==3&&has_field(p.at("compiled_source_sha256"),"examples/n83_rho_matched_eval.rs")&&has_field(p.at("compiled_source_sha256"),"src/cryptanalysis/koblitz_strong_rho.rs")&&has_field(p.at("compiled_source_sha256"),"src/cryptanalysis/semaev_decomp.rs"),"Rho compiled-source binding set mismatch");for(const auto&e:p.at("compiled_source_sha256").object){match(cwd/e.first,e.second.scalar);++source_checks;}require(p.at("compiled_source_sha256").at("examples/n83_rho_matched_eval.rs").scalar==p.at("worker_source_sha256").scalar,"Rho worker source metadata mismatch");
    }
    if(!value_present(raw,"recovered_scalar")){++nonanswers;rows.push_back(row);continue;}
    row.answer=true;row.scalar=raw.at("recovered_scalar").scalar;U d=num(raw.at("recovered_scalar"));require(d>0&&d<instance.r&&mul(d,instance.g)==f->point,"Independent canonical scalar replay failed");row.replay=true;require(uint64_t(d)==f->d,"After-answer validation sidecar mismatch");row.sidecar=true;
    require(truth(raw.at("group_verified"))&&raw.at("exit_code").scalar=="0","Producer group verification/exit flag false");row.online=finite_nonnegative(raw,"target_ms");
    require(truth(observed.at("verified"))&&observed.at("status").scalar=="COMPLETED_VERIFIED"&&std::abs(finite_nonnegative(observed,"online_wall_s")-row.online/1000.)<1e-12,"Observation target verified status or online cost mismatch");
    if(r.arm=="ic"){
     const auto&indices=raw.at("point_indices").array;require(indices.size()==4&&raw.at("x_codes").array.size()==4&&raw.at("pinned_intermediates").array.size()==2,"IC lacks complete four-point witness");P sum;std::vector<P> terms;
     for(size_t k=0;k<4;++k){size_t at=size_t(num(indices[k]));require(at<instance.base.size(),"IC witness outside K600");auto p=instance.base[at];require(p.x==num(raw.at("x_codes").at(k)),"IC witness x-code mismatch");terms.push_back(p);sum=add(sum,p);}
     require(sum==f->point,"Independent IC four-point sum failed");auto left=add(terms[0],terms[1]),right=add(terms[2],terms[3]);require(!left.inf&&!right.inf&&left.x==num(raw.at("pinned_intermediates").at(0))&&right.x==num(raw.at("pinned_intermediates").at(1)),"Independent IC intermediate replay failed");row.witness=true;++witnesses;
     for(const auto&key:{"target_query_ms","target_pdp_ms","target_relation_check_ms","target_descent_ms","target_recovery_check_ms"}){row.phase_sum+=finite_nonnegative(raw,key);}
     require(std::abs(row.phase_sum-finite_nonnegative(raw,"target_phase_sum_ms"))<1e-6,"IC recorded five-phase sum mismatch");
    }else{
     require(raw.at("status").scalar=="PASS_VERIFIED_PUBLIC_TARGET"&&truth(raw.at("verified"))&&truth(raw.at("target_on_curve"))&&truth(raw.at("target_in_subgroup")),"Rho success/membership flags inconsistent");
     U ns=0;for(const auto&key:{"target_query_and_validation","walk_and_collision","generic_scalar_replay"})ns+=num(raw.at("online_phase_ns").at(key));require(ns==num(raw.at("online_phase_sum_ns"))&&ns+num(raw.at("online_timer_boundary_residual_ns"))==num(raw.at("online_interval_ns")),"Rho phase/timer residual mismatch");row.phase_sum=double(ns)/1e6;
    }
    row.gap=row.online-row.phase_sum;row.complete=process_ok&&metadata_error.empty();if(r.arm=="ic")++ic_answers;else ++rho_answers;
    if(r.kind=="aa")require(aa_records.emplace(r.id,raw).second,"Duplicate A/A target record");
   }catch(const std::exception&e){row.error=e.what();++invalid;}rows.push_back(row);
  }
  if(semantic.size()>r.targets.size()){invalid+=semantic.size()-r.targets.size();errors.push_back(r.id+": unpaired extra semantic rows");}
 }
 bool aa_semantic_gate=true;size_t aa_semantic_pairs=0;
 for(size_t pair=1;pair<=5;++pair){
  auto a=aa_records.find(run_id("AA",pair,"A")),b=aa_records.find(run_id("AA",pair,"B"));
  if(a==aa_records.end()||b==aa_records.end()){aa_semantic_gate=false;errors.push_back("Missing verified A/A pair "+std::to_string(pair));continue;}
  bool same=true;for(const auto&key:{"published_q","target","walk_seed","jump_seed","optimization_config","recovered_scalar","outcome","preparation_charges_for_this_seed_only"})same&=equal_json(a->second.at(key),b->second.at(key));
  if(same)++aa_semantic_pairs;else{aa_semantic_gate=false;++invalid;errors.push_back("A/A recorded semantic/counter mismatch pair "+std::to_string(pair));}
 }
 for(const auto&e:files::directory_iterator(campaign/"runs"))if(e.is_directory()&&!planned_dirs.count(files::weakly_canonical(e.path()))){++invalid;errors.push_back("Unplanned run directory: "+e.path().string());}
 for(const auto&e:files::recursive_directory_iterator(campaign))if(e.is_regular_file())hash(e.path());
 auto isolation=campaign.parent_path()/"isolation";if(files::is_directory(isolation))for(const auto&e:files::recursive_directory_iterator(isolation))if(e.is_regular_file())hash(e.path());
 hash(protocol_file);hash(manifest_file);hash(fresh/"generator_receipt.json");hash(fresh/"validation_sidecar.json");
 bool complete=processes==32&&completed==32&&checked==218&&missing==0&&invalid==0&&nonanswers==0&&ic_answers==104&&rho_answers==114&&witnesses==104&&aa_semantic_gate&&aa_semantic_pairs==5&&errors.empty();
 std::ostringstream out;out<<std::setprecision(17)<<"{\n  \"kind\":\"native_n83_matched_answer_and_file_verification\",\n  \"status\":"<<qstr(complete?"PASS":(invalid?"FAIL_INVALID_EVIDENCE":"INCOMPLETE_OR_RETAINED_FAILURES"))<<",\n";
 out<<"  \"campaign_manifest_sha256\":"<<qstr(hash(manifest_file))<<",\n  \"protocol_sha256\":"<<qstr(hash(protocol_file))<<",\n  \"frozen_protocol_sha256\":"<<qstr(hash(protocol_file))<<",\n  \"fixture_manifest_sha256\":"<<qstr(hash(fresh/"generator_receipt.json"))<<",\n  \"validation_sidecar_sha256\":"<<qstr(hash(fresh/"validation_sidecar.json"))<<",\n";
 out<<"  \"curve_slug\":\"icv1-f2m83-t6151469093347-9eba6cd5\",\n  \"planned_processes\":32,\n  \"observed_processes\":"<<processes<<",\n  \"completed_verified_processes\":"<<completed<<",\n  \"planned_target_rows\":218,\n  \"checked_target_rows\":"<<checked<<",\n  \"independent_scalar_replays\":"<<ic_answers+rho_answers<<",\n  \"IC_scalar_replays\":"<<ic_answers<<",\n  \"rho_scalar_replays\":"<<rho_answers<<",\n  \"IC_four_point_witness_replays\":"<<witnesses<<",\n  \"missing_target_rows\":"<<missing<<",\n  \"nonanswer_rows\":"<<nonanswers<<",\n  \"invalid_rows\":"<<invalid<<",\n";
 out<<"  \"recorded_hash_bindings_checked\":"<<hash_checks<<",\n  \"distinct_files_with_recorded_binding\":"<<bound.size()<<",\n  \"all_files_snapshotted\":"<<snapshots.size()<<",\n  \"source_bindings_checked\":"<<source_checks<<",\n  \"clean_environment_records_checked\":"<<env_checks<<",\n  \"resource_envelopes_checked\":"<<resource_checks<<",\n  \"observation_execution_bindings_checked\":"<<observation_checks<<",\n";
 out<<"  \"initial_excluded_signed_orbits\":617,\n  \"fresh_fixtures_reproduced\":32,\n  \"both_prior_eight_fixture_streams_reproduced\":true,\n  \"retained_base_points_and_labels_verified\":99600,\n  \"checked_aa_semantic_pairs\":"<<aa_semantic_pairs<<",\n  \"aa_semantic_gate\":"<<(aa_semantic_gate?"true":"false")<<",\n  \"aa_semantic_scope\":\"Same point, seed, policy, recovered scalar, recorded outcome and preparation charges; not a full intermediate-state trace\",\n  \"complete_panel_gate\":"<<(complete?"true":"false")<<",\n";
 out<<"  \"verification_cost_boundary\":\"This external C++ audit is outside all solver process and online intervals; its resource receipt is reported separately. Both methods already charge built-in group verification.\",\n  \"checker_source_sha256\":"<<qstr(digest(__FILE__))<<",\n  \"checker_executable_sha256\":"<<qstr(digest(argv[0]))<<",\n";
 out<<"  \"not_established\":[\"Independent timing reproduction\",\"Independent host or physical-host exclusivity\",\"Full IC rank and factor-log transcript\",\"Full compiled dependency closure\",\"Common calibrated operation units\",\"Exact five-phase IC online conformance: residual and input-validation mismatch retained\"],\n  \"gate_errors\":[";
 for(size_t i=0;i<errors.size();++i){if(i)out<<",";out<<qstr(errors[i]);}out<<"],\n  \"rows\":[\n";
 for(size_t i=0;i<rows.size();++i){const auto&r=rows[i];if(i)out<<",\n";out<<"    {\"run_id\":"<<qstr(r.run)<<",\"fixture_id\":"<<qstr(r.id)<<",\"arm\":"<<qstr(r.arm)<<",\"process_status\":"<<qstr(r.status)<<",\"row_present\":"<<(r.present?"true":"false")<<",\"answer_present\":"<<(r.answer?"true":"false")<<",\"independent_scalar_replay\":"<<(r.replay?"true":"false")<<",\"IC_four_point_witness_replay\":"<<(r.witness?"true":"false")<<",\"sidecar_checked_after_answer\":"<<(r.sidecar?"true":"false")<<",\"complete_verified_row\":"<<(r.complete?"true":"false")<<",\"recovered_scalar\":"<<(r.scalar.empty()?"null":qstr(r.scalar))<<",\"online_ms\":";if(r.replay)out<<r.online;else out<<"null";out<<",\"phase_sum_ms\":";if(r.replay)out<<r.phase_sum;else out<<"null";out<<",\"phase_gap_ms\":";if(r.replay)out<<r.gap;else out<<"null";out<<",\"error\":"<<(r.error.empty()?"null":qstr(r.error))<<"}";}
 out<<"\n  ],\n  \"files\":[\n";size_t count=0;for(const auto&e:snapshots){if(count++)out<<",\n";out<<"    {\"path\":"<<qstr(e.first.string())<<",\"sha256\":"<<qstr(e.second)<<",\"matched_recorded_binding\":"<<(bound.count(e.first)?"true":"false")<<"}";}out<<"\n  ]\n}\n";create_new(output,out.str());
 std::cout<<(complete?"PASS":"INCOMPLETE/INVALID")<<": processes="<<processes<<"; target_rows="<<checked<<"; independent scalar replays="<<ic_answers+rho_answers<<"; IC witnesses="<<witnesses<<"; invalid="<<invalid<<"; recorded hashes="<<hash_checks<<"\n";return complete?0:(invalid?1:2);
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}

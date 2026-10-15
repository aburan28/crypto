// Read-only diagnosis of v2's lexical JSON binding failures. No solver runs.
#define main frozen_v2_main
#include "../../verify_matched_campaign.cpp"
#undef main

struct Difference {std::string run,comparison,path,left,right;bool numeric=false,equal_double=false,within=false;double ulps=0;};
std::vector<Difference> differences;
void compare(const J&a,const J&b,const std::string&run,const std::string&where,const std::string&path) {
 if(a.type!=b.type||a.array.size()!=b.array.size()||a.object.size()!=b.object.size()){
  differences.push_back({run,where,path,a.type,b.type});return;
 }
 if(a.scalar!=b.scalar){
  Difference d{run,where,path,a.scalar,b.scalar};
  bool floats=a.type=="literal"&&a.scalar.find_first_of(".eE")!=std::string::npos&&b.scalar.find_first_of(".eE")!=std::string::npos;
  if(floats){try{double x=real(a),y=real(b);d.numeric=std::isfinite(x)&&std::isfinite(y);if(d.numeric){d.equal_double=x==y;d.within=same_serialized_f64(x,y);double m=std::max(std::abs(x),std::abs(y));double step=std::nextafter(m,INFINITY)-m;d.ulps=step>0?std::abs(x-y)/step:0;}}catch(const std::exception&){} }
  differences.push_back(d);
 }
 for(size_t i=0;i<a.array.size();++i)compare(a.array[i],b.array[i],run,where,path+"/"+std::to_string(i));
 for(const auto&e:a.object){auto p=b.object.find(e.first);if(p==b.object.end()){differences.push_back({run,where,path+"/"+e.first,"present","missing"});continue;}compare(e.second,p->second,run,where,path+"/"+e.first);}
}
int main(int argc,char**argv){try{
 require(argc==3,"usage: diagnose_numeric_bindings CAMPAIGN NEW_RECEIPT");files::path campaign=argv[1],output=argv[2];require(!files::exists(output),"Refusing existing diagnosis");
 size_t runs=0;
 for(const auto&r:independent_matrix()){
  auto dir=campaign/"runs"/r.id;auto observed=read_json((dir/"observation.json").string()),execution=read_json((dir/"execution.json").string());++runs;
  compare(observed.at("raw_process_receipt"),execution,r.id,"embedded_execution","");
  for(const auto&key:{"elapsed_process_wall_s","cpu_user_s","cpu_system_s","cpu_total_s","peak_rss_bytes"})compare(observed.at(key),execution.at(key),r.id,"process_cost",std::string("/")+key);
  auto rows=lines(dir/"targets.jsonl");for(size_t i=0;i<rows.size();++i)compare(observed.at("target_rows").at(i).at("rawtargetrecord"),rows[i],r.id,"embedded_target",std::string("/")+std::to_string(i));
  auto stdout_rows=lines(dir/"stdout.jsonl");for(const auto&s:stdout_rows)if(has_field(s,"kind")&&s.at("kind").scalar==(r.arm=="rho"?"n83_rho_matched_summary":"compact_orbit_dlp_summary"))compare(observed.at("summary"),s,r.id,"embedded_summary","");
 }
 size_t numeric=0,exact=0,within=0;double max_ulps=0;std::map<std::string,size_t> paths;
 for(const auto&d:differences){numeric+=d.numeric;exact+=d.equal_double;within+=d.within;max_ulps=std::max(max_ulps,d.ulps);++paths[d.path];}
 std::ostringstream out;out<<std::setprecision(17)<<"{\n\"kind\":\"native_readonly_json_binding_diagnosis\",\"checked_processes\":"<<runs<<",\"lexical_mismatches\":"<<differences.size()<<",\"both_finite_f64_tokens\":"<<numeric<<",\"same_decoded_f64\":"<<exact<<",\"within_existing_8epsilon_bound\":"<<within<<",\"maximum_decoded_ulp_distance\":"<<max_ulps<<",\"differences\":[\n";
 for(size_t i=0;i<differences.size();++i){const auto&d=differences[i];if(i)out<<",\n";out<<"{\"run_id\":"<<qstr(d.run)<<",\"comparison\":"<<qstr(d.comparison)<<",\"path\":"<<qstr(d.path)<<",\"left_token\":"<<qstr(d.left)<<",\"right_token\":"<<qstr(d.right)<<",\"both_finite_f64_tokens\":"<<(d.numeric?"true":"false")<<",\"same_decoded_f64\":"<<(d.equal_double?"true":"false")<<",\"within_existing_8epsilon_bound\":"<<(d.within?"true":"false")<<",\"decoded_ulp_distance\":"<<d.ulps<<"}";}
 out<<"\n],\"source_sha256\":"<<qstr(digest(__FILE__))<<",\"executable_sha256\":"<<qstr(digest(argv[0]))<<"}\n";create_new(output,out.str());
 std::cout<<"processes="<<runs<<" lexical_mismatches="<<differences.size()<<" f64="<<numeric<<" exact_double="<<exact<<" within_8eps="<<within<<" max_ulp="<<max_ulps<<"\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}

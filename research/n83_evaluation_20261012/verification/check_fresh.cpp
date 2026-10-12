// Native arithmetic replay only. Run after all timed panel children finish.
// The archived implementation is reused for field/curve arithmetic and JSON
// parsing; its old main is never invoked and archived receipts are untouched.
#define main archived_checker_entry
#include "check.cpp"
#undef main
#include <filesystem>

namespace nfs = std::filesystem;
bool has(const J&j,const std::string&k){return j.object.count(k)!=0;}
bool present(const J&j,const std::string&k){return has(j,k)&&j.at(k).scalar!="null";}
std::string quote(const std::string&s){std::string z="\"";for(char c:s){if(c=='"'||c=='\\')z+='\\';if(c=='\n')z+="\\n";else if(c=='\r')z+="\\r";else if(c=='\t')z+="\\t";else z+=c;}return z+'"';}
std::vector<J> json_lines(const nfs::path&path){std::ifstream f(path);if(!f)throw std::runtime_error("Cannot read "+path.string());std::vector<J> rows;std::string line;while(std::getline(f,line)){if(line.find_first_not_of(" \t\r")==std::string::npos)continue;Parser p;p.s=line;auto j=p.value();p.ws();require(p.i==p.s.size(),"Trailing JSONL data");rows.push_back(std::move(j));}return rows;}
std::vector<P> verify_base(const J&base,uint64_t r,uint64_t lambda){
 check_modulus(base);require(num(base.at("n"))==83&&num(base.at("a"))==1&&num(base.at("subgroup_order"))==r,"Base curve mismatch");
 const auto&raw=base.at("factor_base_point_coordinates").array;const auto&labels=base.at("factor_base_point_labels").array;const auto&reps=base.at("factor_base_representatives").array;
 require(raw.size()==99600&&labels.size()==99600&&reps.size()==600,"Base size mismatch");std::vector<P> points;points.reserve(raw.size());std::set<std::pair<U,U>> distinct;
 for(size_t c=0;c<600;++c){P rep=point(reps[c]);require(on(rep)&&!rep.inf&&mul(r,rep).inf,"Representative subgroup failure");P current=rep;uint64_t coefficient=1;for(int k=0;k<83;++k){for(int sign=0;sign<2;++sign){size_t at=c*166+k*2+sign;P stored=point(raw[at]);require(on(stored)&&stored==(sign?neg(current):current),"Orbit point mismatch");require(num(labels[at].at(0))==c&&num(labels[at].at(1))==(sign?r-coefficient:coefficient),"Orbit coefficient mismatch");require(distinct.insert({stored.x,stored.y}).second,"Duplicate base point");points.push_back(stored);}current=tau(current);coefficient=mm(coefficient,lambda,r);}require(current==rep&&coefficient==1,"Orbit closure failure");}
 return points;
}
struct ReplayRow{std::string file,kind,status,error,recovered,execution_status;bool same_target=false,scalar_replay=false,sidecar_match=false,witness_replay=false;bool answer_present=false;};
int main(int argc,char**argv){try{
 require(argc==4,"usage: check_fresh ARCHIVED_FIXTURE_BASE_DIR FRESH_PANEL_DIR NEW_RECEIPT_JSON");
 const nfs::path archived=argv[1],fresh=argv[2],output=argv[3];require(nfs::exists(fresh),"Fresh panel directory does not exist");require(!nfs::exists(output),"Refusing to overwrite an existing receipt");
 const auto fixture=read_json((archived/"frozen/fixture.json").string());check_modulus(fixture);require(num(fixture.at("n"))==83&&num(fixture.at("a"))==1,"Wrong fixture");
 U x=2;for(int i=0;i<83;++i)x=fs(x);require(x==2&&pgcd(fs(2)^2,MOD)==1,"Field irreducibility failure");uint64_t r=uint64_t(num(fixture.at("subgroup_order")));require(prime(r),"r not prime");
 __int128 t0=2,t1=1;for(int n=2;n<=83;++n){auto t=t1-2*t0;t0=t1;t1=t;}require(TOP+1-U(t1)==num(fixture.at("curve_order"))&&num(fixture.at("cofactor"))*r==num(fixture.at("curve_order")),"Curve recurrence/order mismatch");
 P g=point(fixture.at("generator")),q=point(fixture.at("public_target"));require(on(g)&&on(q)&&mul(r,g).inf&&mul(r,q).inf,"Fixture subgroup mismatch");
 const auto oldrho=read_json((archived/"runs/N83A1K600W12We202610061R21R1/rho.jsonl").string());uint64_t lambda=uint64_t(num(oldrho.at("lambda")));require(lambda!=1&&power(lambda,83,r)==1&&mul(lambda,g)==tau(g),"Eigenvalue mismatch");
 const auto base=read_json((archived/"base_n83_K600.jsonl").string());const auto points=verify_base(base,r,lambda);
 std::vector<nfs::path> producers,executions;for(const auto&entry:nfs::recursive_directory_iterator(fresh)){if(!entry.is_regular_file())continue;auto name=entry.path().filename().string();if(name=="producer.json"||name=="target.jsonl")producers.push_back(entry.path());else if(name=="execution.json")executions.push_back(entry.path());}std::sort(producers.begin(),producers.end());std::sort(executions.begin(),executions.end());
 std::vector<ReplayRow> results;size_t invalid=0,answers=0,rho_rows=0,ic_rows=0;
 for(const auto&path:producers){std::vector<J> rows;try{rows=path.filename()=="target.jsonl"?json_lines(path):std::vector<J>{read_json(path.string())};}catch(const std::exception&e){ReplayRow rr;rr.file=path.string();rr.kind="unparsed_producer_record";rr.status="INVALID_JSON";rr.error=e.what();++invalid;results.push_back(rr);continue;}
  for(const auto&row:rows){ReplayRow rr;rr.file=path.string();const bool ic=path.filename()=="target.jsonl";rr.kind=ic?"IC_target":"strong_rho_producer";if(ic)++ic_rows;else ++rho_rows;rr.status=has(row,"status")?row.at("status").scalar:(present(row,"recovered_scalar")?"ANSWER_RECORDED":"NO_ANSWER_RECORDED");
   auto execution=path.parent_path()/"execution.json";if(nfs::exists(execution)){try{auto e=read_json(execution.string());if(has(e,"status"))rr.execution_status=e.at("status").scalar;}catch(const std::exception&e){rr.error=e.what();}}
   try{
    // Setup/input failures are preserved; they have no target or answer to certify.
    if(!has(row,"generator")||!has(row,"published_q")){require(!present(row,"recovered_scalar"),"Answer lacks exact generator/public-target fields");results.push_back(rr);continue;}
    require(point(row.at("generator"))==g&&point(row.at("published_q"))==q,"Exact generator/public target mismatch");require(num(row.at("n"))==83&&num(row.at("a"))==1,"Curve shape mismatch");rr.same_target=true;
    if(ic){if(has(row,"target"))require(point(row.at("target"))==q,"IC target mismatch");if(has(row,"published_fixture_scalar"))require(row.at("published_fixture_scalar").scalar=="null","IC known scalar input recorded");}
    else{check_modulus(row);require(num(row.at("subgroup_order"))==r&&num(row.at("lambda"))==lambda,"Strong subgroup/eigenvalue mismatch");if(has(row,"public_target_q"))require(point(row.at("public_target_q"))==q,"Strong public target mismatch");if(has(row,"known_answer_input"))require(!truth(row.at("known_answer_input")),"Strong known answer input");if(has(row,"validation_sidecar_opened"))require(!truth(row.at("validation_sidecar_opened")),"Strong producer opened validation sidecar");}
    if(!present(row,"recovered_scalar")){results.push_back(rr);continue;}
    rr.answer_present=true;rr.recovered=row.at("recovered_scalar").scalar;U scalar=num(row.at("recovered_scalar"));require(mul(scalar,g)==q,"Independent scalar replay failed");rr.scalar_replay=true;
    // Read/compare the fixture's validation-only scalar only after answer replay.
    U sidecar=num(fixture.at("fixture_scalar_validation_only"));require(scalar%r==sidecar%r,"Validation-only sidecar answer mismatch");rr.sidecar_match=true;
    if(ic){const auto&idx=row.at("point_indices").array;require(idx.size()==4,"IC answer lacks four-point witness");P total;for(size_t i=0;i<4;++i){size_t at=size_t(num(idx[i]));require(at<points.size(),"IC witness index out of range");if(has(row,"x_codes"))require(points[at].x==num(row.at("x_codes").at(i)),"IC witness x-code mismatch");total=add(total,points[at]);}require(total==q,"IC four-point witness does not sum to Q");if(has(row,"pinned_intermediates")){P a=add(points[size_t(num(idx[0]))],points[size_t(num(idx[1]))]);P b=add(points[size_t(num(idx[2]))],points[size_t(num(idx[3]))]);require(a.x==num(row.at("pinned_intermediates").at(0))&&b.x==num(row.at("pinned_intermediates").at(1)),"IC pinned intermediates mismatch");}rr.witness_replay=true;require(has(row,"group_verified")&&truth(row.at("group_verified")),"IC producer replay flag false");}
    else{require(has(row,"generic_scalar_replay_verified")&&truth(row.at("generic_scalar_replay_verified")),"Strong producer replay flag false");require(rr.status=="PASS_VERIFIED_PUBLIC_TARGET","Strong answer status not verified");}
    ++answers;
   }catch(const std::exception&e){rr.error=e.what();++invalid;}
   results.push_back(rr);
  }
 }
 std::ostringstream disposition;size_t process_count=0;for(const auto&path:executions){auto row=read_json(path.string());if(process_count++)disposition<<",\n";disposition<<"    {\"file\":"<<quote(path.string())<<",\"status\":"<<quote(has(row,"status")?row.at("status").scalar:"MISSING_STATUS")<<",\"return_code\":"<<(has(row,"return_code")?row.at("return_code").scalar:"null")<<",\"termination_signal\":"<<(has(row,"termination_signal")?row.at("termination_signal").scalar:"null")<<"}";}
 std::ofstream out(output);require(bool(out),"Could not create fresh receipt");out<<"{\n  \"kind\":\"native_independent_fresh_panel_answer_replay\",\n  \"status\":"<<quote(invalid?"FAIL":(producers.empty()?"NO_FRESH_PRODUCER_ROWS":"PASS"))<<",\n  \"scope\":\"All retained producer outcomes preserved; standalone arithmetic replay after timed execution, no timing claim\",\n  \"archived_fixture_base_directory\":"<<quote(archived.string())<<",\n  \"fresh_panel_directory\":"<<quote(fresh.string())<<",\n  \"field_curve_subgroup_revalidated\":true,\n  \"factor_base_representatives_verified\":600,\n  \"factor_base_points_and_labels_verified\":99600,\n  \"strong_producer_rows\":"<<rho_rows<<",\n  \"IC_target_rows\":"<<ic_rows<<",\n  \"independently_verified_answers\":"<<answers<<",\n  \"invalid_rows\":"<<invalid<<",\n  \"rows\":[\n";
 for(size_t i=0;i<results.size();++i){const auto&row=results[i];if(i)out<<",\n";out<<"    {\"file\":"<<quote(row.file)<<",\"arm\":"<<quote(row.kind)<<",\"producer_status\":"<<quote(row.status)<<",\"process_status\":"<<quote(row.execution_status)<<",\"answer_present\":"<<(row.answer_present?"true":"false")<<",\"recovered_scalar\":"<<(row.recovered.empty()?"null":quote(row.recovered))<<",\"same_exact_G_Q\":"<<(row.same_target?"true":"false")<<",\"independent_scalar_replay\":"<<(row.scalar_replay?"true":"false")<<",\"sidecar_validation_after_answer\":"<<(row.sidecar_match?"true":"false")<<",\"IC_four_point_witness_replay\":"<<(row.witness_replay?"true":"false")<<",\"error\":"<<(row.error.empty()?"null":quote(row.error))<<"}";}
 out<<"\n  ],\n  \"process_dispositions\":[\n"<<disposition.str()<<"\n  ],\n  \"unverified\":[\"Full factor-log database and 600 rank-equation history\",\"Independent reproduction of timings or operation counts\",\"Independent host/isolation promotion gates\"]\n}\n";out.close();std::cout<<"Fresh archive replay: producer files="<<producers.size()<<"; answers="<<answers<<"; invalid="<<invalid<<"; process dispositions="<<process_count<<"\n";return invalid?1:(producers.empty()?2:0);
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}

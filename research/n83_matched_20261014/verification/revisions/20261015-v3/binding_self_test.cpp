#define main verifier_v3_main
#include "verify_matched_campaign.cpp"
#undef main

int main(int argc,char**argv){try{
 require(argc==2||argc==3,"usage: binding_self_test NEW_RECEIPT [CAMPAIGN]");files::path output=argv[1];require(!files::exists(output),"Refusing existing numeric self-test receipt");
 size_t positive=0,negative=0,actual=0;
 auto check=[&](bool expected,const std::string&a,const std::string&b,BindingContext c,const std::string&label){bool observed=equal_reserialized_json(json_text(a),json_text(b),c);require(observed==expected,("Numeric binding self-test failed: "+label).c_str());if(expected)++positive;else ++negative;};
 auto number=[](double x){std::ostringstream s;s<<std::setprecision(17)<<x;std::string out=s.str();if(out.find_first_of(".eE")==std::string::npos)out+=".0";return out;};
 double base=2100000.,one=std::nextafter(base,INFINITY),two=std::nextafter(one,INFINITY);
 auto target=[](const std::string&n){return std::string("{\"target_ms\":")+n+"}";};
 check(true,R"({"semantic":{"producer_phase_gap_ms":0.00019999999949504857}})",R"({"semantic":{"producer_phase_gap_ms":0.00019999999949504854}})",BindingContext::Process,"actual one-ULP residual");
 check(true,R"({"timing_ms":{"targets_total":3723.018252}})",R"({"timing_ms":{"targets_total":3723.0182520000003}})",BindingContext::Summary,"actual one-ULP summary");
 check(true,target("1e0"),target("1.0"),BindingContext::Target,"different spellings of same f64");
 check(true,target(number(base)),target(number(one)),BindingContext::Target,"one ULP allowed");
 check(false,target(number(base)),target(number(two)),BindingContext::Target,"two ULP rejected");
 check(false,target(number(base)),target(number(base+1e-6)),BindingContext::Target,"one nanosecond at cap in milliseconds rejected");
 check(false,R"({"elapsed_process_wall_s":2100.0})",std::string("{\"elapsed_process_wall_s\":")+number(2100.+1e-9)+"}",BindingContext::Process,"one nanosecond at cap in seconds rejected");
 check(false,R"({"online_interval_ns":9007199254740992})",R"({"online_interval_ns":9007199254740993})",BindingContext::Target,"one integer above2pow53 rejected");
 check(false,R"({"online_phase_ns":{"walk_and_collision":340282366920938463463374607431768211454}})",R"({"online_phase_ns":{"walk_and_collision":340282366920938463463374607431768211455}})",BindingContext::Target,"u128 phase integer one-unit change rejected");
 check(false,R"({"outcome":{"walk_steps":9007199254740992}})",R"({"outcome":{"walk_steps":9007199254740993}})",BindingContext::Target,"counter remains exact");
 check(false,target("1"),target("1.0"),BindingContext::Target,"integer-to-float rejected");
 check(false,target("\"1.0\""),target("1.0"),BindingContext::Target,"string-to-float rejected");
 check(false,R"({"unexpected":{"target_ms":1.0}})",R"({"unexpected":{"target_ms":1.0000000000000002}})",BindingContext::Target,"same leaf under wrong ancestor rejected");
 check(false,R"({"recovered_scalar":"9007199254740992"})",R"({"recovered_scalar":"9007199254740993"})",BindingContext::Target,"scalar string remains exact");
 check(false,R"({"binary_sha256":"abc"})",R"({"binary_sha256":"abd"})",BindingContext::Process,"hash remains exact");
 check(false,R"({"walk_seed":202610146501})",R"({"walk_seed":202610146502})",BindingContext::Target,"seed remains exact");
 check(false,R"({"verified":true})",R"({"verified":false})",BindingContext::Target,"boolean remains exact");
 check(false,R"({"target_ms":1.0})",R"({"target_ms":1.0,"extra":0})",BindingContext::Target,"extra key rejected");
 check(false,R"({"target_ms":1.0})",R"({})",BindingContext::Target,"missing key rejected");
 check(false,R"({"point_indices":[1,2,3,4]})",R"({"point_indices":[2,1,3,4]})",BindingContext::Target,"witness array reorder rejected");
 check(false,R"({"target_ms":1.0})",R"({"target_ms":1.0junk})",BindingContext::Target,"numeric prefix plus junk rejected");
 check(false,R"({"target_ms":1.0})",R"({"target_ms":1e9999})",BindingContext::Target,"overflow rejected");
 check(false,R"({"target_ms":1.0})",R"({"target_ms":NaN})",BindingContext::Target,"NaN rejected");
 for(const auto&bad:{"01.0","+1.0","1.",".1","1e","1e+","-","true","null","Infinity"}){require(!json_float_token(bad),"Malformed token was admitted");++negative;}
 require(equal_json(json_text(R"({"x":1.0})"),json_text(R"({"x":1.0000000000000002})"))==false,"Original exact comparator changed");++negative;
 if(argc==3){files::path campaign=argv[2];for(const auto&r:independent_matrix()){
  auto dir=campaign/"runs"/r.id;auto o=read_json((dir/"observation.json").string()),e=read_json((dir/"execution.json").string());
  require(equal_reserialized_json(o.at("raw_process_receipt"),e,BindingContext::Process),("Actual process binding still mismatches: "+r.id).c_str());++actual;
  for(const auto&key:{"elapsed_process_wall_s","cpu_user_s","cpu_system_s","cpu_total_s","peak_rss_bytes"}){require(equal_reserialized_json(o.at(key),e.at(key),BindingContext::Process,{key}),"Actual cost binding mismatches");++actual;}
  auto rows=lines(dir/"targets.jsonl");for(size_t i=0;i<rows.size();++i){require(equal_reserialized_json(o.at("target_rows").at(i).at("rawtargetrecord"),rows[i],BindingContext::Target),("Actual target binding still mismatches: "+r.id).c_str());++actual;}
  for(const auto&s:lines(dir/"stdout.jsonl"))if(has_field(s,"kind")&&s.at("kind").scalar==(r.arm=="rho"?"n83_rho_matched_summary":"compact_orbit_dlp_summary")){require(equal_reserialized_json(o.at("summary"),s,BindingContext::Summary),("Actual summary binding still mismatches: "+r.id).c_str());++actual;}
 }}
 std::ostringstream out;out<<"{\"kind\":\"native_context_scoped_one_ulp_binding_self_test\",\"status\":\"PASS\",\"positive_cases\":"<<positive<<",\"negative_cases\":"<<negative<<",\"actual_copy_bindings_checked\":"<<actual<<",\"two_ulp_rejected\":true,\"one_nanosecond_rejected_in_seconds_and_milliseconds\":true,\"integer_u128_scalar_seed_hash_bool_key_array_bindings_exact\":true,\"old_exact_comparator_unchanged\":true,\"source_sha256\":"<<qstr(digest(__FILE__))<<",\"numeric_helper_sha256\":"<<qstr(digest(files::path(__FILE__).parent_path()/"numeric_binding.hpp"))<<",\"executable_sha256\":"<<qstr(digest(argv[0]))<<"}\n";create_new(output,out.str());std::cout<<"PASS: positive="<<positive<<" negative="<<negative<<" actual_bindings="<<actual<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}

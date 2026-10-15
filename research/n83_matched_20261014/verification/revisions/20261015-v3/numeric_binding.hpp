#pragma once
// This comparator is used only at four observation reserialization bindings.
// Protocol/configuration/matrix/AA comparisons still use exact equal_json.
enum class BindingContext { Process, Target, Summary };

inline bool json_float_token(const std::string&s) {
 size_t i=0;if(i<s.size()&&s[i]=='-')++i;if(i==s.size())return false;
 if(s[i]=='0')++i;else {if(s[i]<'1'||s[i]>'9')return false;while(i<s.size()&&s[i]>='0'&&s[i]<='9')++i;}
 bool floating=false;
 if(i<s.size()&&s[i]=='.'){floating=true;++i;size_t begin=i;while(i<s.size()&&s[i]>='0'&&s[i]<='9')++i;if(i==begin)return false;}
 if(i<s.size()&&(s[i]=='e'||s[i]=='E')){floating=true;++i;if(i<s.size()&&(s[i]=='+'||s[i]=='-'))++i;size_t begin=i;while(i<s.size()&&s[i]>='0'&&s[i]<='9')++i;if(i==begin)return false;}
 return floating&&i==s.size();
}
inline bool decode_finite_float(const J&j,double&out) {
 if(j.type!="literal"||!json_float_token(j.scalar))return false;
 try{size_t consumed=0;out=std::stod(j.scalar,&consumed);return consumed==j.scalar.size()&&std::isfinite(out);}catch(const std::exception&){return false;}
}
inline bool in_names(const std::string&value,std::initializer_list<const char*>names) {
 return std::any_of(names.begin(),names.end(),[&](const char*n){return value==n;});
}
inline bool canonical_index(const std::string&s) {
 if(s.empty()||(s.size()>1&&s[0]=='0'))return false;
 return std::all_of(s.begin(),s.end(),[](char c){return c>='0'&&c<='9';});
}
inline bool one_f64_roundtrip_step(double a,double b) {
 return std::isfinite(a)&&std::isfinite(b)&&(a==b||std::nextafter(a,b)==b);
}
inline bool float_metadata_path(BindingContext context,const std::vector<std::string>&p) {
 if(context==BindingContext::Target){
  return (p.size()==1&&in_names(p[0],{"target_ms","online_total_ms","online_ms","target_generation_ms_excluded","target_query_ms","target_pdp_ms","target_relation_check_ms","target_descent_ms","target_recovery_check_ms","target_phase_sum_ms"}))
   ||(p.size()==2&&p[0]=="outcome"&&p[1]=="ideal_steps");
 }
 if(context==BindingContext::Summary){
  return (p.size()==1&&in_names(p[0],{"rank_probes_mean","producer_preparation_ms","setup_ms"}))
   ||(p.size()==2&&p[0]=="timing_ms"&&in_names(p[1],{"base_load","normal_basis_and_selftest","index_build","setup_total","rank_stage","linear_algebra","targets_total","target_median","process_total","shared_curve_normal_basis","all_G_only_jumps","targets_online_sum"}));
 }
 if(p.size()==1&&in_names(p[0],{"elapsed_process_wall_s","cpu_user_s","cpu_system_s","cpu_total_s"}))return true;
 if(p.empty()||p[0]!="semantic")return false;
 std::vector<std::string> q(p.begin()+1,p.end());
 if(q.size()>=2&&q[0]=="summary_validation"&&q[1]=="selected_summary")return float_metadata_path(BindingContext::Summary,{q.begin()+2,q.end()});
 if(q.size()>=1&&q[0]=="selected_row")return float_metadata_path(BindingContext::Target,{q.begin()+1,q.end()});
 if(q.size()>=3&&q[0]=="per_target_rows"&&canonical_index(q[1]))q={q.begin()+2,q.end()};
 if(q.size()>=1&&q[0]=="selected_row")return float_metadata_path(BindingContext::Target,{q.begin()+1,q.end()});
 if(q.size()==1&&in_names(q[0],{"target_ms","producer_online_ms","producer_phase_sum_ms","producer_phase_gap_ms","producer_preparation_ms"}))return true;
 return q.size()==2&&q[0]=="producer_phases_ms"&&in_names(q[1],{"target_query_ms","target_pdp_ms","target_relation_check_ms","target_descent_ms","target_recovery_check_ms","walk_ms","collision_ms","recovery_check_ms","validation_ms","target_query_and_validation","walk_and_collision","generic_scalar_replay"});
}
inline bool equal_reserialized_json(const J&a,const J&b,BindingContext context,std::vector<std::string>path={}) {
 if(a.type!=b.type||a.array.size()!=b.array.size()||a.object.size()!=b.object.size())return false;
 if(a.scalar!=b.scalar){double x=0,y=0;if(!float_metadata_path(context,path)||!decode_finite_float(a,x)||!decode_finite_float(b,y)||!one_f64_roundtrip_step(x,y))return false;}
 for(size_t i=0;i<a.array.size();++i){auto child=path;child.push_back(std::to_string(i));if(!equal_reserialized_json(a.array[i],b.array[i],context,child))return false;}
 for(const auto&e:a.object){auto p=b.object.find(e.first);if(p==b.object.end())return false;auto child=path;child.push_back(e.first);if(!equal_reserialized_json(e.second,p->second,context,child))return false;}
 return true;
}

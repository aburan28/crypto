use serde_json::{json,Value};
use std::time::Instant;
include!("../ffi.rs");
fn controls(path:&str)->Value{
 let begin=Instant::now(); let mut s=ffi::Native::new(path,265,8);
 for j in 0..256usize {let x=1+(j%8)as i32;let y=1+((j/8)%8)as i32;let z=9+j as i32;
 assert!(s.clause(&[-x,-y,z]));assert!(s.clause(&[x,-z]));assert!(s.clause(&[y,-z]));}
 let row=(1..=265u32).collect::<Vec<_>>();assert!(s.xor(&row,false));assert!(s.clause(&[-265]));
 let expected=(0..256usize).filter(|x|x.count_ones()%2==0).collect::<std::collections::BTreeSet<_>>();
 let mut got=std::collections::BTreeSet::new();let mut complete=false;
 while begin.elapsed().as_secs_f64()<20.0 {let(v,m)=s.solve((20.0-begin.elapsed().as_secs_f64()).max(0.001));if v==0{complete=true;break;}assert_eq!(v,1,"native interface returned UNKNOWN/error");
 let mask=(0..8usize).fold(0usize,|x,i|x|((m[1+i]==1)as usize)<<i);
 for j in 0..256usize{assert_eq!(m[9+j]==1,((mask>>(j%8))&1!=0)&&((mask>>((j/8)%8))&1!=0));}
 assert_eq!(m[265],0);let parity=m[1..].iter().filter(|x|**x==1).count()%2;assert_eq!(parity,0);assert!(got.insert(mask));
 let block=(1..=8usize).map(|i|if m[i]==1{-(i as i32)}else{i as i32}).collect::<Vec<_>>();if !s.clause(&block){complete=true;break;}}
 assert!(complete);assert_eq!(got,expected);let stats=s.stats();
 let mut empty=ffi::Native::new(path,1,1);assert!(!empty.xor(&[],true));
 json!({"status":"passed","complete_primary_models":got.len(),"independent_oracle_assignments":256,"long_xor_row_arity":265,"all_declared_coordinates_replayed":true,"empty_true_parity_rejected":true,"native_stats":stats,"scope":"interface acceptance and exact small model-set oracle; not the actual compressed norm circuit or a speed comparison"})
}
fn main(){let a=std::env::args().collect::<Vec<_>>();let start=Instant::now();let e=controls(a.get(2).expect("library"));let id=a.get(1).expect("id");let mut r=json!({"schema":"scientific_result_v1","experiment_id":id,"status":"passed","summary":"Long-XOR functional-primary native interface matched an independent complete small model-set oracle.","evidence":{"controls":e,"original_root_count_exact":null,"release_gate_claim":false,"novelty_claim":false},"next_actions":["Encode the exact affine-cone norm representation with this independently tested interface, verify all tiny norm models and freeze a new bounded hard plan.","Keep original S3 ghost/auxiliary coverage and genuine matched F4/F5/F6 performance gates open."]});let mut out=r.to_string();if start.elapsed().as_secs()>=25{r["status"]=json!("inconclusive");r["evidence"]["controls"]["complete_primary_models"]=Value::Null;r["evidence"]["controls"]["status"]=json!("inconclusive");out=r.to_string();}println!("{out}");}

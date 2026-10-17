#[allow(dead_code)]
mod model {
 include!("../../hybrid_norm_xor_20261017_v1/src/main.rs");
 include!("symmetry.rs");
}
fn main(){println!("{}",model::controls());}

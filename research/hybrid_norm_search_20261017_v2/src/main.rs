#[allow(dead_code)]
mod engine {
 include!("../../hybrid_norm_xor_20261017_v1/src/main.rs");
 include!("search.rs");
}
fn main(){engine::entry();}

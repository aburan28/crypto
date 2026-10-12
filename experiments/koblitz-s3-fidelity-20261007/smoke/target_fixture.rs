use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::binary_ecc::curve::{scalar_mul, BinaryPoint};
use num_bigint::BigUint;
fn main() {
  let c=KoblitzCurve::new(0,41).unwrap();
  let d=BigUint::from(123212651130u64);
  match scalar_mul(&c.curve,c.generator(),&d) {
    BinaryPoint::Affine{x,y}=>println!("[{},{}]",x.to_biguint(),y.to_biguint()),
    BinaryPoint::Infinity=>panic!("zero target"),
  }
}

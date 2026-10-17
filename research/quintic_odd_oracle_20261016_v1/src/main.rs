#[allow(dead_code)]
mod candidate {
 include!("../../quintic_resultant_20261016_v1/src/main.rs");
 pub fn check()->Value {
  let mut tested=0;let mut ghosts=0;
  for f in fields().into_iter().filter(|f|matches!(f.k,Kind::Odd(..))){
   let q=f.order().to_usize().unwrap();let l=f.product(&[n(0),n(1)]);
   let products:Vec<P>=(0..=5).map(|k|f.product(&(0..5).map(|i|if i<k{n(1)}else{n(0)}).collect::<Vec<_>>())).collect();
   for mut v in 0..512usize{
    let mut h=vec![];for _ in 0..5{h.push(N::from(v%q));v/=q;}h.push(n(1));
    let expected=oracle_sf(&f,&h);
    for method in 0..=2{assert_eq!(f.sf(&h,method),expected);assert_eq!(f.support(&h,&l,method),products.contains(&h));}tested+=1;
   }
   for a in 0..q{for b in 0..q{
    let quad=vec![N::from(a),N::from(b),n(1)];
    let irreducible=(0..q).all(|x|{let x=N::from(x);!f.add(&f.add(&f.mul(&x,&x),&f.mul(&N::from(b),&x)),&N::from(a)).is_zero()});
    if irreducible{let h=f.mulp(&f.mulp(&quad,&quad),&vec![n(1),n(1)]);for method in 0..=2{assert!(!f.sf(&h,method));assert!(!f.support(&h,&l,method));}ghosts+=1;}
   }}
  }
  json!({"status":"passed","odd_extension_monic_cases":tested,"repeated_nonrational_quadratic_cases":ghosts,"fields":["GF9","GF25"],"source_not_modified":true})
 }
}
fn main(){println!("{}",candidate::check());}

//! Exact known-answer padding control, compiled against the frozen public modules.
//! For m=X^d, every remainder is precisely the first d coefficients.
#![allow(dead_code)]
#[path="../../../tools/isogeny-cli/src/verification_field.rs"]
mod verification_field;
mod field {
    pub type Field=crate::verification_field::PrimeField<3>;
    pub type Fe=crate::verification_field::Element<3>;
}
#[path="../../p192_p224_large_degree_isogenies_20261009/verification_poly.rs"]
mod poly;
use num_bigint::BigUint;
fn main() {
    let p=BigUint::parse_bytes(b"6277101735386680763835789423207666416083908700390324961279",10).unwrap();
    let f=field::Field::new(&p).unwrap();
    let mut checks=0usize;
    for length in [127usize,128,129,130] {
        let d=length-1;
        let mut modulus=vec![f.zero();length];
        modulus[d]=f.one();
        let context=poly::FixedModulus::new(&f,&modulus).unwrap();
        assert_eq!(context.has_reciprocal(),length>=128);
        for n in [0,length-1,length,length+1,length+2,2*length-1,2*length,3*length+3] {
            let dividend=vec![f.one();n];
            let expected=vec![f.one();n.min(d)];
            assert_eq!(context.rem(&dividend),expected);
            assert_eq!(context.rem(&dividend),poly::rem(&f,&dividend,&modulus));
            assert_eq!(context.rem(&dividend),poly::divrem(&f,&dividend,&modulus).1);
            checks+=1;
        }
    }
    println!("{{\"schema\":\"fixed-modulus-zero-tail-control/v1\",\"status\":\"PASS\",\"field\":\"P-192\",\"modulus_lengths\":[127,128,129,130],\"known_answer_reductions\":{checks},\"reciprocal\":\"1 followed by L-1 zero coefficients\",\"scientific_run\":false,\"performance_claim\":false}}");
}


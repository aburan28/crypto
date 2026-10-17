//! Fixed-input panel compiling the repository's actual binary-field source.

#[path = "../../../src/binary_ecc/f2m.rs"]
#[allow(dead_code)]
mod f2m;

use f2m::{F2mElement, IrreduciblePoly};
use std::hint::black_box;
use std::time::Instant;

fn next(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn input(state: &mut u64, m: u32) -> F2mElement {
    let words: Vec<u64> = (0..m.div_ceil(64)).map(|_| next(state)).collect();
    F2mElement::from_words(&words, m)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let iterations: usize = args.get(1).map_or(Ok(200_000), |s| s.parse()).unwrap();
    let rounds: usize = args.get(2).map_or(Ok(5), |s| s.parse()).unwrap();
    assert!(iterations > 0 && rounds > 0);
    let fields = [
        ("f2m113", IrreduciblePoly::deg_113()),
        ("f2m131", IrreduciblePoly::deg_131()),
        ("f2m163", IrreduciblePoly::deg_163()),
        ("f2m233", IrreduciblePoly::deg_233()),
    ];
    println!("field,round,iterations,nanoseconds,checksum");
    for (name, irr) in fields {
        let mut seed = 0x0123_4567_89ab_cdefu64;
        let a: Vec<_> = (0..256).map(|_| input(&mut seed, irr.degree)).collect();
        let b: Vec<_> = (0..256).map(|_| input(&mut seed, irr.degree)).collect();
        for i in 0..256 {
            assert_eq!(a[i].mul(&b[i], &irr), a[i].schoolbook_mul(&b[i], &irr));
        }
        for round in 0..rounds {
            let mut checksum = 0u64;
            let start = Instant::now();
            for i in 0..iterations {
                let j = i & 255;
                let product = black_box(a[j].mul(black_box(&b[j]), black_box(&irr)));
                checksum = checksum.rotate_left(1) ^ product.raw_bits()[0];
            }
            let elapsed = start.elapsed().as_nanos();
            println!("{name},{round},{iterations},{elapsed},{checksum:016x}");
        }
    }
}

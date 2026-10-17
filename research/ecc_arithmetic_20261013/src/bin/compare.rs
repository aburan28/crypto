//! Alternating, same-process comparison of frozen and current field sources.

#[path = "../../baseline_src/f2m.rs"]
#[allow(dead_code)]
mod baseline;
#[path = "../../../../src/binary_ecc/f2m.rs"]
#[allow(dead_code)]
mod candidate;

use std::hint::black_box;
use std::time::Instant;

fn next(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn input(state: &mut u64, m: u32) -> Vec<u64> {
    (0..m.div_ceil(64)).map(|_| next(state)).collect()
}

fn bench_baseline(
    a: &[baseline::F2mElement],
    b: &[baseline::F2mElement],
    irr: &baseline::IrreduciblePoly,
    iterations: usize,
) -> (u128, u64) {
    let mut checksum = 0u64;
    let start = Instant::now();
    for i in 0..iterations {
        let j = i & 255;
        let product = black_box(a[j].mul(black_box(&b[j]), black_box(irr)));
        checksum = checksum.rotate_left(1) ^ product.raw_bits()[0];
    }
    (start.elapsed().as_nanos(), checksum)
}

fn bench_candidate(
    a: &[candidate::F2mElement],
    b: &[candidate::F2mElement],
    irr: &candidate::IrreduciblePoly,
    iterations: usize,
) -> (u128, u64) {
    let mut checksum = 0u64;
    let start = Instant::now();
    for i in 0..iterations {
        let j = i & 255;
        let product = black_box(a[j].mul(black_box(&b[j]), black_box(irr)));
        checksum = checksum.rotate_left(1) ^ product.raw_bits()[0];
    }
    (start.elapsed().as_nanos(), checksum)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let iterations: usize = args.get(1).map_or(Ok(500_000), |s| s.parse()).unwrap();
    let rounds: usize = args.get(2).map_or(Ok(24), |s| s.parse()).unwrap();
    assert!(iterations > 0 && rounds > 0);
    println!("field,round,order,arm,iterations,nanoseconds,checksum");
    for name in ["f2m113", "f2m131", "f2m163", "f2m233"] {
        let (m, base_irr, cand_irr) = match name {
            "f2m113" => (
                113,
                baseline::IrreduciblePoly::deg_113(),
                candidate::IrreduciblePoly::deg_113(),
            ),
            "f2m131" => (
                131,
                baseline::IrreduciblePoly::deg_131(),
                candidate::IrreduciblePoly::deg_131(),
            ),
            "f2m163" => (
                163,
                baseline::IrreduciblePoly::deg_163(),
                candidate::IrreduciblePoly::deg_163(),
            ),
            _ => (
                233,
                baseline::IrreduciblePoly::deg_233(),
                candidate::IrreduciblePoly::deg_233(),
            ),
        };
        let mut seed = 0x0123_4567_89ab_cdefu64;
        let a_words: Vec<_> = (0..256).map(|_| input(&mut seed, m)).collect();
        let b_words: Vec<_> = (0..256).map(|_| input(&mut seed, m)).collect();
        let base_a: Vec<_> = a_words
            .iter()
            .map(|w| baseline::F2mElement::from_words(w, m))
            .collect();
        let base_b: Vec<_> = b_words
            .iter()
            .map(|w| baseline::F2mElement::from_words(w, m))
            .collect();
        let cand_a: Vec<_> = a_words
            .iter()
            .map(|w| candidate::F2mElement::from_words(w, m))
            .collect();
        let cand_b: Vec<_> = b_words
            .iter()
            .map(|w| candidate::F2mElement::from_words(w, m))
            .collect();
        for i in 0..256 {
            let want = base_a[i].schoolbook_mul(&base_b[i], &base_irr);
            assert_eq!(base_a[i].mul(&base_b[i], &base_irr), want);
            assert_eq!(
                cand_a[i].mul(&cand_b[i], &cand_irr).raw_bits(),
                want.raw_bits()
            );
        }
        for round in 0..rounds {
            let (first, second) = if round & 1 == 0 {
                ("baseline", "candidate")
            } else {
                ("candidate", "baseline")
            };
            let order = if round & 1 == 0 { "AB" } else { "BA" };
            for arm in [first, second] {
                let (elapsed, checksum) = if arm == "baseline" {
                    bench_baseline(&base_a, &base_b, &base_irr, iterations)
                } else {
                    bench_candidate(&cand_a, &cand_b, &cand_irr, iterations)
                };
                println!("{name},{round},{order},{arm},{iterations},{elapsed},{checksum:016x}");
            }
        }
    }
}

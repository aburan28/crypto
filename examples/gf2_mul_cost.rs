//! The instruction cost of one `Gf2` multiplication and one squaring, as the
//! strong rho calls them: from ordinary (not `pclmulqdq`-compiled) code, so
//! each product is a call to the folded multiply.
//!
//! `gf2_mul_cost n muls sqrs` runs `muls` dependent multiplications, then
//! `sqrs` dependent squarings, in the field the strong-rho fixtures use at
//! `n` (53 or 61), and prints a checksum.  Two Callgrind runs that differ
//! only in `muls` (or `sqrs`) give the instructions per operation as the
//! difference over the extra count: everything else in the process is the
//! same.
use crypto_lib::binary_ecc::IrreduciblePoly;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use std::hint::black_box;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: u32 = args.get(1).map_or(53, |s| s.parse().expect("n"));
    let muls: u64 = args.get(2).map_or(1_000_000, |s| s.parse().expect("muls"));
    let sqrs: u64 = args.get(3).map_or(1_000_000, |s| s.parse().expect("sqrs"));
    let low_terms = match n {
        53 => vec![0, 1, 2, 6],
        61 => vec![0, 1, 2, 5],
        _ => panic!("n must be 53 or 61"),
    };
    let gf = Gf2::new(&IrreduciblePoly {
        degree: n,
        low_terms,
    });
    let mask = (1u64 << n) - 1;
    let x = black_box(0x0123_4567_89ab_cdef & mask);
    let mut acc = black_box(0x0fed_cba9_8765_4321 & mask);
    for _ in 0..muls {
        acc = gf.mul(acc, x);
    }
    let mut sq = black_box(acc | 1);
    for _ in 0..sqrs {
        sq = gf.sqr(sq);
    }
    println!(
        "{{\"n\":{n},\"muls\":{muls},\"sqrs\":{sqrs},\"checksum\":\"{:016x}\"}}",
        acc ^ sq
    );
}

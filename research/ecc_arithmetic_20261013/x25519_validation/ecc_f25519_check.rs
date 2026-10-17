#[path = "../../../src/ecc/f25519.rs"]
mod f25519;

use f25519::Fe;

fn hex(bytes: &[u8; 32]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn emit(a: [u8; 32], b: [u8; 32]) {
    let fa = Fe::from_bytes(&a);
    let fb = Fe::from_bytes(&b);
    println!(
        "{} {} {} {} {} {}",
        hex(&a),
        hex(&b),
        hex(&fa.add(&fb).to_bytes()),
        hex(&fa.sub(&fb).to_bytes()),
        hex(&fa.mul(&fb).to_bytes()),
        hex(&fa.invert().to_bytes())
    );
}

fn main() {
    let zero = [0u8; 32];
    let mut one = zero;
    one[0] = 1;
    let mut p = [0xffu8; 32];
    p[0] = 0xed;
    p[31] = 0x7f;
    let mut p_minus_one = p;
    p_minus_one[0] = 0xec;
    let mut p_plus_one = p;
    p_plus_one[0] = 0xee;
    let mut high_bit = zero;
    high_bit[31] = 0x80;
    for (a, b) in [
        (zero, zero),
        (one, zero),
        (p_minus_one, one),
        (p, one),
        (p_plus_one, p_minus_one),
        ([0xff; 32], p),
        (high_bit, one),
        ([0xff; 32], [0xff; 32]),
    ] {
        emit(a, b);
    }

    let mut seed = 0x91e1_0da5_c79e_7b1du64;
    for _ in 0..992 {
        let mut a = [0u8; 32];
        let mut b = [0u8; 32];
        for bytes in [&mut a, &mut b] {
            for chunk in bytes.chunks_exact_mut(8) {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                chunk.copy_from_slice(&seed.to_le_bytes());
            }
        }
        emit(a, b);
    }
}

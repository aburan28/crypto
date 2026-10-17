#[path = "../../../src/ecc/f25519.rs"]
mod f25519;

use f25519::Fe;

fn hex(bytes: &[u8; 32]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn main() {
    let mut seed = 0x91e1_0da5_c79e_7b1du64;
    for _ in 0..1000 {
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
}

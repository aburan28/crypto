#[path = "../../../src/ecc/p521_field.rs"]
mod p521_field;
use p521_field::P521Fe;
fn h(bytes: &[u8]) -> String { bytes.iter().map(|b| format!("{b:02x}")).collect() }
fn main() {
    let mut state = 0x41b7_6d8e_579d_1531u64;
    for i in 0..1000 {
        let mut a = [0u8;66]; let mut b = [0u8;66];
        for out in [&mut a, &mut b] {
            for chunk in out.chunks_mut(8) {
                state ^= state << 13; state ^= state >> 7; state ^= state << 17;
                chunk.copy_from_slice(&state.to_le_bytes()[..chunk.len()]);
            }
        }
        if i == 0 { a = [0;66]; b = [0;66]; }
        if i == 1 { a = [255;66]; b = [255;66]; }
        if i == 2 { a = [255;66]; a[65] = 1; b[0] = 1; }
        let x=P521Fe::from_le_bytes(&a); let y=P521Fe::from_le_bytes(&b);
        println!("{} {} {} {} {} {} {} {}",h(&a),h(&b),h(&x.to_le_bytes()),h(&y.to_le_bytes()),h(&x.add(&y).to_le_bytes()),h(&x.sub(&y).to_le_bytes()),h(&x.mul(&y).to_le_bytes()),h(&x.invert().to_le_bytes()));
    }
}

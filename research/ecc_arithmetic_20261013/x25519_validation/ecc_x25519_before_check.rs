#[path = "ecc_before.rs"]
mod ecc;

fn hex(bytes: &[u8; 32]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn hex32(value: &str) -> [u8; 32] {
    let mut bytes = [0u8; 32];
    for (i, byte) in bytes.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&value[i * 2..i * 2 + 2], 16).unwrap();
    }
    bytes
}

fn main() {
    let vectors = [
        (
            "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4",
            "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c",
            "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552",
        ),
        (
            "4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d",
            "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493",
            "95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957",
        ),
    ];
    for (scalar, u, expected) in vectors {
        assert_eq!(
            ecc::x25519::x25519(&hex32(scalar), &hex32(u)),
            hex32(expected)
        );
    }

    let mut seed = 0x8ba1_267c_341e_d0f5u64;
    for _ in 0..100 {
        let mut scalar = [0u8; 32];
        let mut u = [0u8; 32];
        for bytes in [&mut scalar, &mut u] {
            for chunk in bytes.chunks_exact_mut(8) {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                chunk.copy_from_slice(&seed.to_le_bytes());
            }
        }
        println!("{} {} {}", hex(&scalar), hex(&u), hex(&ecc::x25519::x25519(&scalar, &u)));
    }
}

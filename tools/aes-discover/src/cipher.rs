//! Independent portable AES reference. Bytes use FIPS column-major order.
//! Table lookups and exposed keys are intended for synthetic experiments only.
use serde::{Deserialize, Serialize};
use std::sync::LazyLock;

pub type Block = [u8; 16];

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub key_bits: usize,
    pub rounds: usize,
    pub final_mix_columns: bool,
}

impl Config {
    pub fn check(self) -> Result<(), String> {
        let full = match self.key_bits {
            128 => 10,
            192 => 12,
            256 => 14,
            _ => return Err("key_bits must be 128, 192 or 256".into()),
        };
        if !(1..=full).contains(&self.rounds) {
            return Err(format!(
                "rounds must be in 1..={full} for AES-{}",
                self.key_bits
            ));
        }
        Ok(())
    }
}

pub const fn mul(mut a: u8, mut b: u8) -> u8 {
    let mut p = 0;
    let mut i = 0;
    while i < 8 {
        if b & 1 != 0 {
            p ^= a;
        }
        let high = a >> 7;
        a <<= 1;
        if high != 0 {
            a ^= 0x1b;
        }
        b >>= 1;
        i += 1;
    }
    p
}

pub const fn inverse(x: u8) -> u8 {
    let mut p = 1;
    let mut a = x;
    let mut n = 254;
    while n > 0 {
        if n & 1 != 0 {
            p = mul(p, a);
        }
        a = mul(a, a);
        n >>= 1;
    }
    p // Includes the AES convention inverse(0) = 0.
}

pub static SBOX: LazyLock<[u8; 256]> = LazyLock::new(|| {
    std::array::from_fn(|x| {
        let y = inverse(x as u8);
        y ^ y.rotate_left(1) ^ y.rotate_left(2) ^ y.rotate_left(3) ^ y.rotate_left(4) ^ 0x63
    })
});
pub static INV_SBOX: LazyLock<[u8; 256]> = LazyLock::new(|| {
    let mut out = [0; 256];
    for i in 0..256 {
        out[SBOX[i] as usize] = i as u8;
    }
    out
});

pub fn xor(a: Block, b: Block) -> Block {
    std::array::from_fn(|i| a[i] ^ b[i])
}

pub fn shift(s: &mut Block, backward: bool) {
    let old = *s;
    for c in 0..4 {
        for r in 0..4 {
            let source = if backward {
                (c + 4 - r) % 4
            } else {
                (c + r) % 4
            };
            s[4 * c + r] = old[4 * source + r];
        }
    }
}

pub fn mix(s: &mut Block, backward: bool) {
    for col in s.chunks_exact_mut(4) {
        let a: [u8; 4] = col.try_into().expect("four-byte column");
        if backward {
            for r in 0..4 {
                col[r] = mul(14, a[r])
                    ^ mul(11, a[(r + 1) % 4])
                    ^ mul(13, a[(r + 2) % 4])
                    ^ mul(9, a[(r + 3) % 4]);
            }
        } else {
            let total = a[0] ^ a[1] ^ a[2] ^ a[3];
            for r in 0..4 {
                col[r] = a[r] ^ total ^ mul(2, a[r] ^ a[(r + 1) % 4]);
            }
        }
    }
}

pub fn linear(s: &mut Block, use_mix: bool) {
    shift(s, false);
    if use_mix {
        mix(s, false);
    }
}

#[derive(Clone)]
pub struct Aes {
    pub config: Config,
    pub round_keys: Vec<Block>,
}

impl Aes {
    pub fn new(key: &[u8], config: Config) -> Result<Self, String> {
        config.check()?;
        if key.len() * 8 != config.key_bits {
            return Err("key length differs from config".into());
        }
        let nk = key.len() / 4;
        let full_rounds = nk + 6;
        let mut words: Vec<[u8; 4]> = key.chunks_exact(4).map(|c| c.try_into().unwrap()).collect();
        let mut rcon = 1;
        for i in nk..4 * (full_rounds + 1) {
            let mut t = words[i - 1];
            if i % nk == 0 {
                t.rotate_left(1);
                for x in &mut t {
                    *x = SBOX[*x as usize];
                }
                t[0] ^= rcon;
                rcon = mul(2, rcon);
            } else if nk == 8 && i % nk == 4 {
                for x in &mut t {
                    *x = SBOX[*x as usize];
                }
            }
            words.push(std::array::from_fn(|j| words[i - nk][j] ^ t[j]));
        }
        let round_keys = words
            .chunks_exact(4)
            .map(|w| std::array::from_fn(|i| w[i / 4][i % 4]))
            .collect();
        Ok(Self { config, round_keys })
    }

    pub fn encrypt(&self, input: Block) -> Block {
        let mut s = xor(input, self.round_keys[0]);
        for r in 1..=self.config.rounds {
            for b in &mut s {
                *b = SBOX[*b as usize];
            }
            linear(
                &mut s,
                r < self.config.rounds || self.config.final_mix_columns,
            );
            s = xor(s, self.round_keys[r]);
        }
        s
    }

    pub fn decrypt(&self, input: Block) -> Block {
        let mut s = input;
        for r in (1..=self.config.rounds).rev() {
            s = xor(s, self.round_keys[r]);
            if r < self.config.rounds || self.config.final_mix_columns {
                mix(&mut s, true);
            }
            shift(&mut s, true);
            for b in &mut s {
                *b = INV_SBOX[*b as usize];
            }
        }
        xor(s, self.round_keys[0])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{block, hex, unhex};

    #[test]
    fn fips_128_192_256() {
        let pt = block("00112233445566778899aabbccddeeff").unwrap();
        for (bits, rounds, expected) in [
            (128, 10, "69c4e0d86a7b0430d8cdb78070b4c55a"),
            (192, 12, "dda97ca4864cdfe06eaf70a0ec0d7191"),
            (256, 14, "8ea2b7ca516745bfeafc49904b496089"),
        ] {
            let key: Vec<u8> = (0..bits / 8).map(|i| i as u8).collect();
            let aes = Aes::new(
                &key,
                Config {
                    key_bits: bits,
                    rounds,
                    final_mix_columns: false,
                },
            )
            .unwrap();
            assert_eq!(hex(&aes.encrypt(pt)), expected);
            assert_eq!(aes.decrypt(block(expected).unwrap()), pt);
        }
    }

    #[test]
    fn round_trips_and_linear_inverse() {
        for bits in [128, 192, 256] {
            for rounds in 1..=bits / 32 + 6 {
                for final_mix_columns in [false, true] {
                    let aes = Aes::new(
                        &vec![0x52; bits / 8],
                        Config {
                            key_bits: bits,
                            rounds,
                            final_mix_columns,
                        },
                    )
                    .unwrap();
                    let pt = std::array::from_fn(|i| (17 * i) as u8);
                    assert_eq!(aes.decrypt(aes.encrypt(pt)), pt);
                }
            }
        }
        for x in 1..=255 {
            assert_eq!(mul(x, inverse(x)), 1);
        }
        assert_eq!(SBOX[0], 0x63);
        assert_eq!(SBOX[0x53], 0xed);
        assert!(unhex("éé").is_err());
    }

    #[test]
    fn column_xor_is_preserved_and_two_byte_cancellation_exists() {
        for x in 0..=255 {
            let mut s = [x; 16];
            s[0] ^= 0x71;
            let before = s;
            mix(&mut s, false);
            for c in 0..4 {
                assert_eq!(
                    before[4 * c..4 * c + 4].iter().fold(0, |a, b| a ^ b),
                    s[4 * c..4 * c + 4].iter().fold(0, |a, b| a ^ b)
                );
            }
            mix(&mut s, true);
            assert_eq!(s, before);
        }
        let mut s = [0; 16];
        s[0] = 3;
        s[1] = 2;
        mix(&mut s, false);
        assert_eq!(s[0], 0); // Branch number 5 does not mean every output is active.
    }
}

//! Exact component enumeration and descriptive mechanism checks. Component
//! identities and diffusion plots are not complete-cipher distinguishers.
use crate::cipher::{inverse, mul, xor, Aes, Block, Config, INV_SBOX, SBOX};
use crate::search::{public_key, RandomPermutation};
use crate::{hex, stream};
use rand::seq::SliceRandom;
use rand::RngCore;
use serde::Serialize;

pub fn ddt(sbox: &[u8; 256]) -> Vec<Vec<u16>> {
    let mut t = vec![vec![0u16; 256]; 256];
    for a in 0..256 {
        for x in 0..256 {
            t[a][(sbox[x] ^ sbox[x ^ a]) as usize] += 1;
        }
    }
    t
}

#[derive(Debug, Serialize)]
pub struct TableReport {
    pub schema: String,
    pub ddt_max_nonzero_input: u16,
    pub ddt_row_sums_are_256: bool,
    pub ddt_entries_are_even: bool,
    pub lat_max_abs_walsh_nonzero_masks: i32,
    pub bct_max_nonzero_differences_including_diagonal: u16,
    pub bct_zero_row_and_column_are_256: bool,
    pub enumerated_ddt_cells: usize,
    pub enumerated_lat_cells: usize,
    pub enumerated_bct_cells: usize,
    pub diagonal_bct_max: u16,
}

pub fn tables() -> TableReport {
    let table = ddt(&SBOX);
    let mut lat_max = 0i32;
    for a in 1..256usize {
        for b in 1..256usize {
            let mut walsh = 0i32;
            for x in 0..256usize {
                walsh += if ((a & x).count_ones() ^ (b & SBOX[x] as usize).count_ones()) & 1 == 0 {
                    1
                } else {
                    -1
                };
            }
            lat_max = lat_max.max(walsh.abs());
        }
    }
    let mut bct_max = 0;
    let mut diagonal = 0;
    let mut bct_zero_ok = true;
    for a in 0..256usize {
        for b in 0..256usize {
            let count = (0..256usize)
                .filter(|&x| {
                    INV_SBOX[(SBOX[x] as usize) ^ b] ^ INV_SBOX[(SBOX[x ^ a] as usize) ^ b]
                        == a as u8
                })
                .count() as u16;
            if a != 0 && b != 0 {
                bct_max = bct_max.max(count);
            } else {
                bct_zero_ok &= count == 256;
            }
            if a == b && a != 0 {
                diagonal = diagonal.max(count);
            }
        }
    }
    TableReport {
        schema: "aes-discover/sbox-tables/v1".into(),
        ddt_max_nonzero_input: *table[1..].iter().flatten().max().unwrap(),
        ddt_row_sums_are_256: table
            .iter()
            .all(|r| r.iter().map(|&x| x as u32).sum::<u32>() == 256),
        ddt_entries_are_even: table.iter().flatten().all(|n| n % 2 == 0),
        lat_max_abs_walsh_nonzero_masks: lat_max,
        bct_max_nonzero_differences_including_diagonal: bct_max,
        bct_zero_row_and_column_are_256: bct_zero_ok,
        enumerated_ddt_cells: 65536,
        enumerated_lat_cells: 255 * 255,
        enumerated_bct_cells: 65536,
        diagonal_bct_max: diagonal,
    }
}

#[derive(Serialize)]
pub struct IntegralRow {
    pub key_index: usize,
    pub structure_index: usize,
    pub active_plaintext_byte: usize,
    pub ciphertext_xor: String,
    pub balanced_output_bytes: usize,
    pub permutation_xor: String,
    pub permutation_balanced_output_bytes: usize,
    pub one_member_mutant_xor: String,
}

pub fn integral(
    config: Config,
    seed: u64,
    keys: usize,
    structures: usize,
    active: usize,
) -> Result<serde_json::Value, String> {
    config.check()?;
    if keys == 0 || keys > 1024 || structures == 0 || structures > 4096 || active >= 16 {
        return Err("keys 1..1024, structures 1..4096, active byte 0..15 required".into());
    }
    let mut rows = Vec::new();
    for k in 0..keys {
        let aes = Aes::new(
            &public_key(seed, "integral/keys", k as u64, config.key_bits),
            config,
        )?;
        let mut prp = RandomPermutation::new(seed, "integral/prp", k as u64);
        for j in 0..structures {
            let mut base = [0; 16];
            stream(seed, "integral/bases", (k * structures + j) as u64).fill_bytes(&mut base);
            let mut sum = [0; 16];
            let mut null = [0; 16];
            let mut mutant = [0; 16];
            for v in 0..256 {
                let mut p = base;
                p[active] = v as u8;
                let y = aes.encrypt(p);
                sum = xor(sum, y);
                null = xor(null, prp.encrypt(p));
                if v == 255 {
                    p[(active + 1) % 16] ^= 1;
                }
                mutant = xor(mutant, aes.encrypt(p));
            }
            rows.push(IntegralRow {
                key_index: k,
                structure_index: j,
                active_plaintext_byte: active,
                ciphertext_xor: hex(&sum),
                balanced_output_bytes: sum.iter().filter(|&&x| x == 0).count(),
                permutation_xor: hex(&null),
                permutation_balanced_output_bytes: null.iter().filter(|&&x| x == 0).count(),
                one_member_mutant_xor: hex(&mutant),
            });
        }
    }
    Ok(
        serde_json::json!({"schema":"aes-discover/integral/v1","cipher":config,"seed":seed,
        "keys":keys,"structures_per_key":structures,"texts_per_structure":256,
        "aes_queries_including_mutants":keys*structures*512,"permutation_queries":keys*structures*256,
        "rows":rows,"scope":"descriptive integral reproduction with consistent PRP and one-member sibling mutant; no new integral or universal-key theorem"}),
    )
}

fn affine(x: u8) -> u8 {
    x ^ x.rotate_left(1) ^ x.rotate_left(2) ^ x.rotate_left(3) ^ x.rotate_left(4)
}
fn ratio(x: [u8; 4]) -> u8 {
    mul(
        mul(x[0] ^ x[2], x[1] ^ x[3]),
        inverse(mul(x[0] ^ x[3], x[1] ^ x[2])),
    )
}

pub fn quartets(seed: u64, trials: u64) -> Result<serde_json::Value, String> {
    if trials == 0 || trials > 100_000_000 {
        return Err("trials must be 1..100000000".into());
    }
    let mut rng = stream(seed, "quartets", 0);
    let mut affine_inv = [0; 256];
    for x in 0..256 {
        affine_inv[affine(x as u8) as usize] = x as u8;
    }
    let mut random_box: [u8; 256] = std::array::from_fn(|i| i as u8);
    random_box.shuffle(&mut rng);
    let mut inversion_matches = 0;
    let mut aes_matches = 0;
    let mut undone_matches = 0;
    let mut random_matches = 0;
    let mut poles_skipped = 0;
    let mut accepted = 0;
    for _ in 0..trials {
        let kin = rng.next_u32() as u8;
        let kout = rng.next_u32() as u8;
        let mut x = [0; 4];
        for i in 0..4 {
            loop {
                let v = rng.next_u32() as u8;
                if !x[..i].contains(&v) {
                    x[i] = v;
                    break;
                }
            }
        }
        if x.contains(&kin) {
            poles_skipped += 1;
            continue;
        }
        accepted += 1;
        let input_ratio = ratio(x);
        let inv = x.map(|v| inverse(v ^ kin) ^ kout);
        let aes = x.map(|v| SBOX[(v ^ kin) as usize] ^ kout);
        inversion_matches += u64::from(input_ratio == ratio(inv));
        aes_matches += u64::from(input_ratio == ratio(aes));
        undone_matches += u64::from(input_ratio == ratio(aes.map(|v| affine_inv[v as usize])));
        random_matches += u64::from(input_ratio == ratio(x.map(|v| random_box[v as usize])));
    }
    Ok(
        serde_json::json!({"schema":"aes-discover/cross-ratio-component/v1","seed":seed,"generated_quartets":trials,
        "accepted_nondegenerate_quartets":accepted,"inversion_poles_skipped":poles_skipped,
        "key_wrapped_inversion_matches":inversion_matches,"raw_aes_sbox_matches":aes_matches,
        "aes_output_linear_affine_undone_matches":undone_matches,"random_bijection_matches":random_matches,
        "formula":"CR(a,b,c,d)=((a+c)*(b+d))/((a+d)*(b+c)) in GF(256), distinct inputs, x xor input_key != 0",
        "scope":"known local projective identity and affine-layer ablation only; no reproduction of the complete 2026/1549 attack or whole-round invariant"}),
    )
}

pub fn schedule(seed: u64, keys: usize, key_bit: usize) -> Result<serde_json::Value, String> {
    if keys == 0 || keys > 1024 || key_bit >= 128 {
        return Err("keys 1..1024 and common key bit 0..127 required".into());
    }
    let mut rows = Vec::new();
    for bits in [128, 192, 256] {
        for k in 0..keys {
            let key = public_key(seed, "schedule/keys", k as u64, bits);
            let mut changed = key.clone();
            changed[key_bit / 8] ^= 1 << (key_bit % 8);
            let config = Config {
                key_bits: bits,
                rounds: bits / 32 + 6,
                final_mix_columns: false,
            };
            let a = Aes::new(&key, config)?;
            let b = Aes::new(&changed, config)?;
            let mut null_rng = stream(
                seed,
                "schedule/independent-round-key-differences",
                (bits * keys + k) as u64,
            );
            for (round, (x, y)) in a.round_keys.iter().zip(&b.round_keys).enumerate() {
                let d = xor(*x, *y);
                let mut null_difference = [0; 16];
                null_rng.fill_bytes(&mut null_difference);
                rows.push(serde_json::json!({"key_bits":bits,"key_index":k,"round":round,
                    "key_difference_bit":key_bit,"active_round_key_bytes":d.iter().filter(|&&x|x!=0).count(),
                    "round_key_hamming_weight":d.iter().map(|x|x.count_ones()).sum::<u32>(),
                    "independent_round_key_active_bytes":null_difference.iter().filter(|&&x|x!=0).count(),
                    "independent_round_key_hamming_weight":null_difference.iter().map(|x|x.count_ones()).sum::<u32>()}));
            }
        }
    }
    Ok(
        serde_json::json!({"schema":"aes-discover/key-schedule/v1","seed":seed,"keys_per_variant":keys,
        "rows":rows,"scope":"one-bit related-key diffusion diagnostic, with independent random round-key difference reference; not a related-key attack or single-key weakness claim"}),
    )
}

pub fn boomerang(
    config: Config,
    seed: u64,
    samples: u64,
    alpha: Block,
    delta: Block,
) -> Result<serde_json::Value, String> {
    config.check()?;
    if samples == 0 || samples > 10_000_000 || alpha == [0; 16] || delta == [0; 16] {
        return Err("nonzero frozen differences and samples 1..10000000 required".into());
    }
    let aes = Aes::new(
        &public_key(seed, "boomerang/key", 0, config.key_bits),
        config,
    )?;
    let mut prp = RandomPermutation::new(seed, "boomerang/prp", 0);
    let mut rng = stream(seed, "boomerang/plaintexts", 0);
    let mut hits = 0;
    let mut control_hits = 0;
    let mut trivial = 0;
    let mut control_trivial = 0;
    for _ in 0..samples {
        let mut p = [0; 16];
        rng.fill_bytes(&mut p);
        let q = xor(p, alpha);
        let (a, b) = (aes.encrypt(p), aes.encrypt(q));
        trivial += u64::from(xor(a, b) == delta);
        hits += u64::from(xor(aes.decrypt(xor(a, delta)), aes.decrypt(xor(b, delta))) == alpha);
        let (a, b) = (prp.encrypt(p), prp.encrypt(q));
        control_trivial += u64::from(xor(a, b) == delta);
        control_hits +=
            u64::from(xor(prp.decrypt(xor(a, delta)), prp.decrypt(xor(b, delta))) == alpha);
    }
    Ok(
        serde_json::json!({"schema":"aes-discover/boomerang/v1","cipher":config,"seed":seed,"samples":samples,
        "frozen_alpha":hex(&alpha),"frozen_delta":hex(&delta),"returns":hits,"permutation_returns":control_hits,
        "trivial_swaps":trivial,"permutation_trivial_swaps":control_trivial,"queries_each_arm":4*samples,
        "scope":"descriptive fixed-difference forward/inverse-oracle experiment; observed ciphertext differences never choose delta; no significance or CP-only claim"}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_aes_component_constants_including_bct_diagonal() {
        let t = tables();
        assert_eq!(t.ddt_max_nonzero_input, 4);
        assert!(t.ddt_row_sums_are_256 && t.ddt_entries_are_even);
        assert_eq!(t.lat_max_abs_walsh_nonzero_masks, 32);
        assert_eq!(t.bct_max_nonzero_differences_including_diagonal, 6);
        assert!(t.bct_zero_row_and_column_are_256);
    }
    #[test]
    fn integral_positive_null_and_nonvacuity_mutant() {
        let r = integral(
            Config {
                key_bits: 128,
                rounds: 3,
                final_mix_columns: true,
            },
            42,
            3,
            2,
            0,
        )
        .unwrap();
        for row in r["rows"].as_array().unwrap() {
            assert_eq!(row["balanced_output_bytes"], 16);
            assert_ne!(row["one_member_mutant_xor"], "00".repeat(16));
            assert_ne!(row["permutation_xor"], "00".repeat(16));
        }
    }
    #[test]
    fn cross_ratio_is_component_scoped() {
        let r = quartets(42, 4096).unwrap();
        assert_eq!(
            r["key_wrapped_inversion_matches"],
            r["accepted_nondegenerate_quartets"]
        );
        assert_eq!(
            r["aes_output_linear_affine_undone_matches"],
            r["accepted_nondegenerate_quartets"]
        );
        assert_ne!(
            r["raw_aes_sbox_matches"],
            r["accepted_nondegenerate_quartets"]
        );
    }
}

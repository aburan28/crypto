//! Deterministic CM fixture construction; certification is delegated to plan().
use isogeny_algos::{curve::{self, Pt}, field::{is_prime, Rng, Zp}};
use std::{env, fs, path::Path};
use volcano_cli::{plan, Input};

fn main() {
    let root_arg = env::args().nth(1).expect("study directory required");
    let root = Path::new(&root_arg);
    for (label, dk, conductor, j, start) in [
        ("two", 7u64, 8192u64, -3375i64, 1u64 << 29),
        ("three", 11, 6561, -32768, 1u64 << 29),
        ("mixed", 19, 1728, -884736, 1u64 << 30),
    ] {
        let mut found = false;
        for offset in 0..100_000u64 {
            let t = start + offset;
            let numerator = t as u128 * t as u128 + dk as u128 * conductor as u128 * conductor as u128;
            if numerator % 4 != 0 { continue; }
            let p = (numerator / 4) as u64;
            if !(1u64 << 51..1u64 << 61).contains(&p) || !is_prime(p) { continue; }
            let field = Zp::new(p);
            let c = curve::from_j(&field, j.rem_euclid(p as i64) as u64);
            let point = curve::random_point_f(&field, &c, &mut Rng::new(1));
            for n in [p + 1 - t, p + 1 + t] {
                if curve::pmul(&field, &c, &point, n as u128) != Pt::Inf { continue; }
                let input = Input { p, a:c.a, b:c.b, order:Some(n) };
                let receipt = match plan(input.clone(), "auto") {
                    Ok(p) => p,
                    Err(e) => { eprintln!("candidate {label} p={p} n={n}: {e}"); continue; }
                };
                assert_eq!(receipt.fundamental_discriminant, -(dk as i64));
                assert_eq!(receipt.frobenius_conductor, conductor);
                fs::write(root.join("inputs").join(format!("{label}.json")), format!("{}\n", serde_json::to_string_pretty(&input).unwrap())).unwrap();
                fs::write(root.join("evidence").join(format!("{label}-plan.json")), format!("{}\n", serde_json::to_string_pretty(&receipt).unwrap())).unwrap();
                println!("{}", serde_json::json!({"case":label,"p":p,"trace":receipt.trace,"conductor":conductor,"vertices":receipt.expected_vertices,"strata":receipt.strata.len(),"identity":receipt.source_identity}));
                found = true;
                break;
            }
            if found { break; }
        }
        assert!(found, "candidate search exhausted for {label}");
    }
}

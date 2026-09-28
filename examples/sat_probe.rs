use crypto_lib::cryptanalysis::sat::{SolveResult, Solver};
use crypto_lib::cryptanalysis::semaev_corpus::CORPUS;
use crypto_lib::cryptanalysis::semaev_sat::{encode_semaev_s4_with, S4Options, XorEncoding};
use std::time::Instant;

struct SplitMix(u64);
impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        ((self.next() as u128 * n as u128) >> 64) as u64
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args[1].as_str() {
        "corpus" => {
            let enc_kind = if args.get(3).map(|s| s.as_str()) == Some("cnf") {
                XorEncoding::Cnf
            } else {
                XorEncoding::Native
            };
            for c in CORPUS.iter().filter(|c| c.name.starts_with(args[2].as_str())) {
                let mut e = encode_semaev_s4_with(
                    c.n,
                    c.l,
                    &c.irr(),
                    &c.b(),
                    &c.x_r(),
                    S4Options {
                        encoding: enc_kind,
                        break_symmetry: true,
                    },
                );
                let t = Instant::now();
                let r = e.solver.solve();
                println!(
                    "{} {:?} {:.1} ms conflicts {}",
                    c.name,
                    r,
                    t.elapsed().as_secs_f64() * 1e3,
                    e.solver.conflicts()
                );
            }
        }
        "r3" => {
            let n: u32 = args[2].parse().unwrap();
            let m: usize = args[3].parse().unwrap();
            for seed in 1..=12u64 {
                let mut rng = SplitMix(seed);
                let mut s = Solver::new(n);
                for _ in 0..m {
                    let mut c: Vec<i32> = Vec::new();
                    while c.len() < 3 {
                        let v = 1 + rng.below(n as u64) as i32;
                        if c.iter().any(|&l: &i32| l.abs() == v) {
                            continue;
                        }
                        c.push(if rng.next() & 1 == 1 { v } else { -v });
                    }
                    s.add_clause(c);
                }
                let t = Instant::now();
                let r = s.solve();
                println!(
                    "seed {seed} {:?} {:.1} ms conflicts {}",
                    r,
                    t.elapsed().as_secs_f64() * 1e3,
                    s.conflicts()
                );
            }
        }
        "xor" => {
            let n: u32 = args[2].parse().unwrap();
            let rows: usize = args[3].parse().unwrap();
            let width: usize = args[4].parse().unwrap();
            let m: usize = args[5].parse().unwrap();
            for seed in 0..6u64 {
                let mut rng = SplitMix(0x5eed_0000 + seed);
                let hidden: Vec<bool> = (0..n).map(|_| rng.next() & 1 == 1).collect();
                let mut s = Solver::new(n);
                for _ in 0..rows {
                    let mut vars: Vec<u32> = Vec::new();
                    while vars.len() < width {
                        let v = 1 + rng.below(n as u64) as u32;
                        if !vars.contains(&v) {
                            vars.push(v);
                        }
                    }
                    let rhs = vars.iter().fold(false, |a, &v| a ^ hidden[(v - 1) as usize]);
                    s.add_xor(&vars, rhs);
                }
                let mut k = 0;
                while k < m {
                    let mut c: Vec<i32> = Vec::new();
                    while c.len() < 3 {
                        let v = 1 + rng.below(n as u64) as i32;
                        if c.iter().any(|&l: &i32| l.abs() == v) {
                            continue;
                        }
                        c.push(if rng.next() & 1 == 1 { v } else { -v });
                    }
                    if c.iter().any(|&l| hidden[(l.abs() - 1) as usize] == (l > 0)) {
                        s.add_clause(c);
                        k += 1;
                    }
                }
                s.conflict_budget = 200_000;
                let t = Instant::now();
                let r = s.solve();
                println!(
                    "seed {seed} {:?} {:.1} ms conflicts {} xor_row_ops {}",
                    r,
                    t.elapsed().as_secs_f64() * 1e3,
                    s.conflicts(),
                    s.stats.xor_row_ops
                );
                let _ = SolveResult::Sat;
            }
        }
        _ => panic!(),
    }
}

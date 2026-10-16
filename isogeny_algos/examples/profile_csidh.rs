//! CSIDH-512 action with sqrt-Velu leaves from various degrees (scratch profiling); also checks
//! that every threshold gives the same codomain.
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::path::csidh::Csidh;
use std::time::Instant;

fn main() {
    let mut cs = Csidh::<FpM<8>>::csidh512();
    let f = cs.fp.clone();
    let mut rng = Rng::new(77);
    let keys: Vec<Vec<i32>> = (0..12)
        .map(|_| (0..74).map(|_| (rng.next() % 11) as i32 - 5).collect())
        .collect();
    let mut reference = vec![];
    for (round, thr) in [
        u64::MAX,
        263,
        151,
        u64::MAX,
        263,
        151,
        u64::MAX,
        263,
        151,
        u64::MAX,
        263,
        151,
    ]
    .iter()
    .enumerate()
    {
        cs.sqrt_velu_min = *thr;
        let mut ts = vec![];
        for (i, key) in keys.iter().enumerate() {
            let mut r = Rng::new(1000 + i as u64);
            let t = Instant::now();
            let a = cs.action_fast(f.zero(), key, &mut r);
            ts.push(t.elapsed().as_secs_f64() * 1e3);
            if round == 0 {
                reference.push(a);
            } else {
                assert_eq!(a, reference[i], "threshold {thr}");
            }
        }
        ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!("threshold {thr:>20}: median {:.1} ms", ts[ts.len() / 2]);
    }
}

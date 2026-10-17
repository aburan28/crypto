//! x^q mod g for a random monic g of degree d over word-size and 127-, 255- and 511-bit prime
//! fields: the table
//! reducer against Newton/Barrett (`poly::POWMOD_TABLE_N` forces the choice), results checked
//! equal. `cargo run --release --example profile_powmod`
use isogeny_algos::field::*;
use isogeny_algos::fpm::FpM;
use isogeny_algos::poly::{self, POWMOD_TABLE_N};
use std::time::Instant;

fn run<F: Field>(name: &str, f: &F) {
    let mut rng = Rng::new(3);
    println!("# {name}: deg  table_ms  barrett_ms  barrett/table");
    for d in [16, 24, 32, 40, 48, 56, 63, 72, 80, 90, 100, 120, 140, 160] {
        let mut g: Vec<F::E> = (0..d).map(|_| f.random(&mut rng)).collect();
        g.push(f.one());
        let x = poly::x_poly(f);
        let mut best = [f64::MAX; 2];
        let mut out = vec![];
        for (k, ov) in [usize::MAX, 0].into_iter().enumerate() {
            POWMOD_TABLE_N.with(|c| c.set(Some(ov)));
            for _ in 0..3 {
                let t = Instant::now();
                let r = poly::powmod_big(f, &x, &f.q(), &g);
                best[k] = best[k].min(t.elapsed().as_secs_f64() * 1e3);
                if out.len() == k {
                    out.push(r);
                }
            }
        }
        POWMOD_TABLE_N.with(|c| c.set(None));
        assert_eq!(out[0], out[1], "{name} deg {d}");
        println!(
            "{d:5} {:9.2} {:11.2} {:8.2}",
            best[0],
            best[1],
            best[1] / best[0]
        );
    }
}

fn main() {
    run("FpM1 (2^61-1)", &FpM::<1>::from_dec("2305843009213693951"));
    run("Zp (40-bit)", &Zp::new(next_prime((1u64 << 39) + 12345)));
    run("Zp (61-bit)", &Zp::new(next_prime((1u64 << 60) + 12345)));
    run(
        "FpM2 (2^127-1)",
        &FpM::<2>::from_dec("170141183460469231731687303715884105727"),
    );
    run(
        "FpM4 (2^255-19)",
        &FpM::<4>::from_dec(
            "57896044618658097711785492504343953926634992332820282019728792003956564819949",
        ),
    );
    let csidh = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";
    run("FpM8 (CSIDH-512)", &FpM::<8>::from_dec(csidh));
}

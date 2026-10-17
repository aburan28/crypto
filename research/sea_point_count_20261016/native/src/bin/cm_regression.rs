//! An untimed regression prompted by the formal sample-only counterexample.
use isogeny_algos::{
    api::PrimeCurve,
    gls_loop::{self, Instance},
    int::Int,
    json::Json,
};

fn exhaustive(c: &PrimeCurve) -> u64 {
    let p = c.p.to_i128().unwrap() as u64;
    let (a, b) = (c.a.mod_u64(p), c.b.mod_u64(p));
    1 + (0..p)
        .map(|x| {
            (0..p)
                .filter(|&y| y * y % p == (x * x % p * x + a * x + b) % p)
                .count() as u64
        })
        .sum::<u64>()
}

fn main() {
    let mut failures = 0;
    for (p, d, t, v, h, r) in [
        (11i64, -7i64, 4i64, 2i64, 4u64, 29i64),
        (11, -7, -4, 2, 4, 29),
        (59, -23, 12, 2, 4, 877),
        (59, -23, -12, 2, 4, 877),
    ] {
        for seed in 1..=8 {
            let inst = Instance {
                d,
                p: Int::from(p),
                t: Int::from(t),
                v: Int::from(v),
                h,
                r: Int::from(r),
                tries: 0,
            };
            let expected = (p + 1 - t) as u64;
            let generated = gls_loop::build_curve(&inst, seed);
            let actual = generated.as_ref().ok().map(|cm| exhaustive(&cm.curve));
            let pass = actual == Some(expected);
            if !pass {
                failures += 1;
            }
            println!(
                "{}",
                Json::obj(vec![
                    ("p", Json::Num(p)),
                    ("d", Json::Num(d)),
                    ("trace", Json::Num(t)),
                    ("seed", Json::Num(seed as i64)),
                    ("expected_order", Json::Num(expected as i64)),
                    (
                        "actual_order",
                        actual.map(|n| Json::Num(n as i64)).unwrap_or(Json::Null)
                    ),
                    ("pass", Json::Bool(pass)),
                    (
                        "error",
                        generated.err().map(Json::str).unwrap_or(Json::Null)
                    )
                ])
                .dump()
            );
        }
    }
    let invalid = Instance {
        d: -23,
        p: Int::from(59i64),
        t: Int::from(-36i64),
        v: Int::from(2i64),
        h: 4,
        r: Int::from(877i64),
        tries: 0,
    };
    let rejected = gls_loop::build_curve(&invalid, 1).is_err();
    if !rejected {
        failures += 1;
    }
    println!(
        "{}",
        Json::obj(vec![
            ("invalid_norm_rejected", Json::Bool(rejected)),
            ("failures", Json::Num(failures))
        ])
        .dump()
    );
    if failures != 0 {
        std::process::exit(1);
    }
}

//! Exact small-eigenvalue-order screening on public standard curve parameters.
use isogeny_algos::{
    api::PrimeCurve,
    bigint::Big,
    field::{is_prime, Zp},
    int::Int,
    json::Json,
};
fn pow(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = r * a % p;
        }
        a = a * a % p;
        e >>= 1;
    }
    r
}
fn order(a: u64, p: u64) -> u64 {
    let mut n = p - 1;
    let mut residual = n;
    let mut d = 2;
    while d * d <= residual {
        if residual % d == 0 {
            while n % d == 0 && pow(a, n / d, p) == 1 {
                n /= d;
            }
            while residual % d == 0 {
                residual /= d;
            }
        }
        d += 1;
    }
    if residual > 1 {
        while n % residual == 0 && pow(a, n / residual, p) == 1 {
            n /= residual;
        }
    }
    n
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 1);
    let mut rows = vec![];
    for (name, n) in [
        (
            "p192",
            "6277101735386680763835789423176059013767194773182842284081",
        ),
        (
            "p224",
            "26959946667150639794667015087019625940457807714424391721682722368061",
        ),
    ] {
        let c = PrimeCurve::preset(name).unwrap();
        let n = Int::from_big(&Big::from_dec(n));
        let t = &(&c.p + &Int::one()) - &n;
        for ell in 1010..=65537u64 {
            if !is_prime(ell) {
                continue;
            }
            let f = Zp::new(ell);
            let tm = t.mod_u64(ell);
            let pm = c.p.mod_u64(ell);
            let disc = (tm * tm + ell - (4 * pm) % ell) % ell;
            if disc == 0 {
                continue;
            }
            let Some(s) = f.sqrt(disc) else {
                continue;
            };
            let roots = [
                (tm + s) % ell * ((ell + 1) / 2) % ell,
                (tm + ell - s) % ell * ((ell + 1) / 2) % ell,
            ];
            let orders = [order(roots[0], ell), order(roots[1], ell)];
            if *orders.iter().min().unwrap() <= 8 {
                rows.push(Json::obj(vec![
                    ("curve", Json::str(name)),
                    ("ell", Json::Num(ell as i64)),
                    (
                        "eigenvalues",
                        Json::Arr(roots.into_iter().map(|v| Json::Num(v as i64)).collect()),
                    ),
                    (
                        "eigenvalue_orders",
                        Json::Arr(orders.into_iter().map(|v| Json::Num(v as i64)).collect()),
                    ),
                ]));
            }
        }
    }
    let out = Json::obj(vec![
        ("schema", Json::str("small-eigenvalue-orders/v1")),
        ("lower", Json::Num(1010)),
        ("upper", Json::Num(65537)),
        ("maximum_selected_order", Json::Num(8)),
        ("candidates", Json::Arr(rows)),
    ]);
    std::fs::write(&args[0], format!("{}\n", out.dump())).unwrap();
    println!("Exact screen written to {}", args[0]);
}

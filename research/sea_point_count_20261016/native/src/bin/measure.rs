use isogeny_algos::{api::{self, PrimeCurve, PrimeFieldInt}, bigint::Big,
    curve::{self, Curve}, field::{Field, Rng, Zp}, find::{sea, modpoly::Phi},
    fpm::FpM, icv1, int::Int, json::Json};
use std::{collections::HashMap, hint::black_box, time::Instant};

fn dec(s: &str) -> Int { Int::from_big(&Big::from_dec(s)) }
fn fixture(case: &str) -> (PrimeCurve, Int, usize) {
    let (p,a,b,n,reps) = match case {
        "bsgs24" => ("10796507","9425520","6283680","10790192", 1000),
        "bsgs32" => ("2790788831","666561205","2453577277","2790782094", 300),
        "bsgs40" | "sea40_cached" => ("775471435903","483575337897","620539964477","775471436947", 100),
        "sea192_cached" | "count192_cold" => {
            return (PrimeCurve::preset("p192").unwrap(), dec("6277101735386680763835789423176059013767194773182842284081"), 1);
        },
        "cm256_control" => {
            return (PrimeCurve::preset("secp256k1").unwrap(),dec("115792089237316195423570985008687907852837564279074904382605163141518161494337"),10);
        },
        _ => panic!("unknown workload: {case}"),
    };
    (PrimeCurve::new(dec(p),dec(a),dec(b)).unwrap(),dec(n),reps)
}

fn timed_sea<F: PrimeFieldInt>(f:&F,c:&PrimeCurve,n:&Int,reps:usize) -> (u128,u128) {
    let e=Curve::new(f.elem(&c.a),f.elem(&c.b));
    let mut cache:HashMap<usize,Phi<F>>=HashMap::new();
    let setup=Instant::now();
    let warm=sea::sea(f,&e,131,&mut cache,&mut Rng::new(1)).expect("warm SEA").0;
    assert_eq!(warm,*n.mag());
    let setup_ns=setup.elapsed().as_nanos();
    let start=Instant::now();
    for i in 0..reps {
        let result=sea::sea(f,&e,131,&mut cache,&mut Rng::new(1+(i%4) as u64)).expect("SEA").0;
        assert_eq!(result,*n.mag());
        black_box(result);
    }
    (start.elapsed().as_nanos(),setup_ns)
}

fn main() {
    let case=std::env::args().nth(1).expect("workload name");
    let (c,n,reps)=fixture(&case);
    let id=icv1::prime(&c.p,&c.a,&c.b,&n).unwrap();
    let (elapsed_ns,setup_ns)=if case.starts_with("bsgs") {
        let f=Zp::new(c.p.to_i128().unwrap() as u64);
        let e=Curve::new(f.elem(&c.a),f.elem(&c.b));
        let start=Instant::now();
        for i in 0..reps {
            let result=curve::order(&f,&e,&mut Rng::new(1+(i%4) as u64));
            assert_eq!(Int::from(result),n);
            black_box(result);
        }
        (start.elapsed().as_nanos(),0)
    } else if case=="sea40_cached" {
        timed_sea(&Zp::new(c.p.to_i128().unwrap() as u64),&c,&n,reps)
    } else if case=="sea192_cached" {
        timed_sea(&FpM::<3>::new(c.p.mag()),&c,&n,reps)
    } else {
        let start=Instant::now();
        for i in 0..reps {
            let result=api::count_points(&c,1+(i%4) as u64).expect("count");
            assert_eq!(result.order,n);
            assert!(result.certificate.certified && result.twist_certificate.certified);
            black_box(result);
        }
        (start.elapsed().as_nanos(),0)
    };
    println!("{}",Json::obj(vec![
        ("schema",Json::str("sea-bench/v1")),("case",Json::str(case)),
        ("curve",Json::str(id.slug)),("order",Json::str(n)),
        ("calls",Json::Num(reps as i64)),("seed_schedule",Json::str("1 + i % 4")),
        ("elapsed_ns",Json::str(elapsed_ns)),("setup_ns",Json::str(setup_ns)),
        ("verified",Json::Bool(true)),
    ]).dump());
}

//! Independent prime-field and curve replay of held-out χ relation solves.

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

type El = [u64; 2];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Point {
    x: El,
    y: El,
    infinity: bool,
}

#[derive(Clone, Copy)]
struct Field {
    p: u64,
    nu: u64,
}

impl Field {
    fn add(self, x: El, y: El) -> El {
        [(x[0] + y[0]) % self.p, (x[1] + y[1]) % self.p]
    }
    fn sub(self, x: El, y: El) -> El {
        [
            (x[0] + self.p - y[0]) % self.p,
            (x[1] + self.p - y[1]) % self.p,
        ]
    }
    fn mul(self, x: El, y: El) -> El {
        let p = self.p as u128;
        let real =
            (x[0] as u128 * y[0] as u128 + self.nu as u128 * x[1] as u128 * y[1] as u128) % p;
        let imag = (x[0] as u128 * y[1] as u128 + x[1] as u128 * y[0] as u128) % p;
        [real as u64, imag as u64]
    }
    fn scale(self, x: El, scalar: u64) -> El {
        self.mul(x, [scalar % self.p, 0])
    }
    fn inv(self, x: El) -> El {
        let norm = self.sub(
            self.mul([x[0], 0], [x[0], 0]),
            self.scale([x[1] * x[1] % self.p, 0], self.nu),
        )[0];
        assert_ne!(norm, 0);
        let inv = pow_mod(norm, self.p - 2, self.p);
        self.mul([x[0], (self.p - x[1]) % self.p], [inv, 0])
    }
    fn neg(self, x: El) -> El {
        self.sub([0, 0], x)
    }
}

fn pow_mod(mut x: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = ((result as u128 * x as u128) % modulus as u128) as u64;
        }
        x = ((x as u128 * x as u128) % modulus as u128) as u64;
        exponent >>= 1;
    }
    result
}

fn identity() -> Point {
    Point {
        x: [0, 0],
        y: [0, 0],
        infinity: true,
    }
}

fn neg(f: Field, point: Point) -> Point {
    Point {
        y: f.neg(point.y),
        ..point
    }
}

fn add(f: Field, a: El, p: Point, q: Point) -> Point {
    if p.infinity {
        return q;
    }
    if q.infinity {
        return p;
    }
    if p.x == q.x && p.y == f.neg(q.y) {
        return identity();
    }
    let slope = if p == q {
        let numerator = f.add(f.scale(f.mul(p.x, p.x), 3), a);
        f.mul(numerator, f.inv(f.scale(p.y, 2)))
    } else {
        f.mul(f.sub(q.y, p.y), f.inv(f.sub(q.x, p.x)))
    };
    let x = f.sub(f.sub(f.mul(slope, slope), p.x), q.x);
    let y = f.sub(f.mul(slope, f.sub(p.x, x)), p.y);
    Point {
        x,
        y,
        infinity: false,
    }
}

fn scalar_mul(f: Field, a: El, point: Point, mut scalar: u64) -> Point {
    let mut result = identity();
    let mut base = point;
    while scalar != 0 {
        if scalar & 1 != 0 {
            result = add(f, a, result, base);
        }
        base = add(f, a, base, base);
        scalar >>= 1;
    }
    result
}

fn on_curve(f: Field, a: El, b: El, p: Point) -> bool {
    p.infinity || f.mul(p.y, p.y) == f.add(f.add(f.mul(f.mul(p.x, p.x), p.x), f.mul(a, p.x)), b)
}

fn point(v: &Value) -> Point {
    if v["infinity"] == true {
        return identity();
    }
    Point {
        x: [v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()],
        y: [v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()],
        infinity: false,
    }
}

fn expected_pairs(f: Field, a: El, signed: &[Point], target: Point) -> BTreeSet<(usize, usize)> {
    let mut out = BTreeSet::new();
    for (i, &p) in signed.iter().enumerate() {
        for (j, &q) in signed.iter().enumerate().skip(i) {
            if add(f, a, p, q) == target {
                out.insert((i, j));
            }
        }
    }
    out
}

fn modmul(a: u64, b: u64) -> u64 {
    a * b % 229
}

fn insert(pivots: &mut [Option<[u64; 10]>; 9], mut row: [u64; 10]) -> bool {
    for col in 0..9 {
        if row[col] == 0 { continue; }
        if let Some(pivot) = pivots[col] {
            let factor = row[col];
            for j in col..10 {
                row[j] = (row[j] + 229 - modmul(factor, pivot[j])) % 229;
            }
        } else {
            let inverse = pow_mod(row[col], 227, 229);
            for item in row.iter_mut().skip(col) { *item = modmul(*item, inverse); }
            pivots[col] = Some(row);
            return true;
        }
    }
    assert_eq!(row[9], 0);
    false
}

fn velu(f: Field, source_a: El, t: Point, p: Point) -> Point {
    if p.infinity || p == t || p == neg(f, t) { return identity(); }
    let plus = add(f, source_a, p, t);
    let minus = add(f, source_a, p, neg(f, t));
    Point {
        x: f.sub(f.add(f.add(p.x, plus.x), minus.x), f.scale(t.x, 2)),
        y: f.add(f.add(p.y, plus.y), minus.y),
        infinity: false,
    }
}

fn check_row(f: Field, public: &Value, checker: &Value, signed: &[Point],
    row: &Value, round: u64, arm: &str, position: usize) -> (u64, u64) {
    let a = [190, 0];
    let g = point(&public["generator"]);
    let mapped_g = point(&public["mapped_generator"]);
    let q = point(&public["targets"][position]["source"]);
    let mapped_q = point(&public["targets"][position]["destination"]);
    assert_eq!(row["schema"], "different_j_chi_heldout_run/v2");
    assert_eq!(row["public_sha256"], checker["public_sha256"]);
    assert_eq!(row["round"], round);
    assert_eq!(row["arm"], arm);
    assert_eq!(row["position"], position);
    assert_eq!(point(&row["target"]), q);
    assert_eq!(point(&row["mapped_target"]), mapped_q);
    let seed = (0x2026_1015u64 ^ 211)
        .wrapping_add(round.wrapping_mul(0x9e37_79b9))
        .wrapping_add(((position + 64) as u64).wrapping_mul(0x0123_4567));
    assert_eq!(row["seed"], seed);
    let scalar = row["recovered"].as_u64().unwrap();
    assert_eq!(scalar, checker["scalars"][position].as_u64().unwrap());
    assert_eq!(scalar_mul(f, a, g, scalar), q);
    assert_eq!(scalar_mul(f, [130,0], mapped_g, scalar), mapped_q);
    assert_eq!(row["verified"], true);
    assert!(row["peak_rss_bytes"].as_u64().unwrap() <= 512*1024*1024);
    let mut matched = 0u64;
    let trials = row["trials"].as_u64().unwrap();
    if arm == "chi-resultant" || arm == "direct-pair" {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut pivots: [Option<[u64; 10]>; 9] = [None; 9];
        let witnesses = row["relation_witnesses"].as_array().unwrap();
        let mut rank = 0u64;
        for trial in 0..trials {
            let ra = rng.gen_range(0..229u64);
            let rb = rng.gen_range(0..229u64);
            let candidate = add(f,a,scalar_mul(f,a,g,ra),scalar_mul(f,a,q,rb));
            let projected = scalar_mul(f,a,candidate,195);
            if projected.infinity { continue; }
            let pairs = expected_pairs(f,a,signed,projected);
            if let Some(&(i,j)) = pairs.first() {
                let witness = &witnesses[matched as usize];
                assert_eq!(witness["a"],ra);
                assert_eq!(witness["b"],rb);
                assert_eq!(witness["i"],i);
                assert_eq!(witness["j"],j);
                assert_eq!(point(&witness["projected"]),projected);
                let mut row_vec = [0u64; 10];
                for index in [i,j] {
                    let col = index/2;
                    row_vec[col] = (row_vec[col] + if index%2==0 {1} else {228}) % 229;
                }
                row_vec[8]=(229-modmul(195,rb))%229;
                row_vec[9]=modmul(195,ra);
                if insert(&mut pivots,row_vec) { rank+=1; }
                matched+=1;
                assert!(rank<9 || trial+1==trials);
            }
        }
        assert_eq!(matched as usize,witnesses.len());
        assert_eq!(row["relations"],matched);
        assert_eq!(row["rank"],rank);
        assert_eq!(rank,9);
        assert!(trials<=1832);
    } else if arm=="source-bsgs" {
        let mut baby=std::collections::BTreeMap::new();
        let mut cursor=identity();
        for small in 0..16u64 {
            baby.insert((cursor.infinity,cursor.x,cursor.y),small);
            cursor=add(f,a,cursor,g);
        }
        let giant_step=neg(f,scalar_mul(f,a,g,16));
        let mut cursor=q;
        let mut found=None;
        for giant in 0..15u64 {
            if let Some(&small)=baby.get(&(cursor.infinity,cursor.x,cursor.y)) {
                found=Some((giant*16+small)%229);
                assert_eq!(trials,giant+1);
                break;
            }
            cursor=add(f,a,cursor,giant_step);
        }
        assert_eq!(found,Some(scalar));
    } else {
        let attempts=row["rho_attempts"].as_array().unwrap();
        assert!(!attempts.is_empty() && attempts.len()<=32);
        assert_eq!(attempts.last().unwrap()["verified"],true);
        for attempt in attempts {
            assert!(attempt["steps"].as_u64().unwrap()<=4580);
            if attempt["verified"]==true {
                assert_eq!(attempt["recovered"],scalar);
            }
        }
    }
    (trials,matched)
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(),4,"usage: checker PUBLIC.json CHECKER.json RUN_DIRECTORY");
    let public_bytes=fs::read(&args[1]).unwrap();
    let checker_bytes=fs::read(&args[2]).unwrap();
    let public:Value=serde_json::from_slice(&public_bytes).unwrap();
    let checker:Value=serde_json::from_slice(&checker_bytes).unwrap();
    assert_eq!(public["schema"],"different_j_chi_heldout_public/v1");
    assert_eq!(checker["schema"],"different_j_chi_heldout_checker/v1");
    assert_eq!(checker["public_sha256"],sha256_hex(&public_bytes));
    assert_eq!(checker["target_count"],4);
    let f=Field{p:211,nu:2};
    let a=[190,0];
    let g=point(&public["generator"]);
    let mapped_g=point(&public["mapped_generator"]);
    let t=Point{x:[3,0],y:[1,0],infinity:false};
    assert!(on_curve(f,a,[37,0],g));
    assert!(on_curve(f,[130,0],[179,0],mapped_g));
    assert_eq!(velu(f,a,t,g),mapped_g);
    let mut signed=Vec::new();
    for (column,row) in public["full_base_witness_source"].as_array().unwrap().iter().enumerate() {
        let raw=point(row);
        assert!(on_curve(f,a,[37,0],raw));
        let projected=scalar_mul(f,a,raw,195);
        assert_eq!(scalar_mul(f,a,projected,229),identity());
        assert_eq!(velu(f,a,t,projected),point(&public["full_base_projected"][column]));
        signed.push(projected);
        signed.push(neg(f,projected));
    }
    let arms=["chi-resultant","direct-pair","source-bsgs","rho-negation","rho-negation-aa"];
    let mut rows=0u64;
    let mut relation_witnesses=0u64;
    let mut trials=0u64;
    for round in 0..3 {
        for arm in arms {
            for position in 0..4 {
                let id=format!("r{round}_{arm}_t{position}");
                let root=Path::new(&args[3]);
                let exit=fs::read_to_string(root.join(format!("{id}.exit"))).unwrap();
                assert_eq!(exit.trim(),"0");
                let raw=fs::read(root.join(format!("{id}.json"))).unwrap();
                let row:Value=serde_json::from_slice(&raw).unwrap();
                let meter:Value=serde_json::from_slice(
                    &fs::read(root.join(format!("{id}.meter.json"))).unwrap()).unwrap();
                assert_eq!(meter["schema"],"different_j_chi_heldout_meter/v1");
                assert_eq!(meter["exit_code"],0);
                assert_eq!(meter["timed_out"],false);
                assert!(meter["child_peak_rss_bytes"].as_u64().unwrap()<=512*1024*1024);
                assert_eq!(meter["stdout_bytes"].as_u64().unwrap(),raw.len() as u64);
                let (this_trials,this_witnesses)=check_row(f,&public,&checker,&signed,
                    &row,round,arm,position);
                trials+=this_trials;
                relation_witnesses+=this_witnesses;
                rows+=1;
            }
        }
    }
    assert_eq!(rows,60);
    println!("{}",json!({"schema":"different_j_chi_heldout_check/v2","status":"pass",
        "rows":rows,"trials":trials,"relation_witnesses":relation_witnesses,
        "public_sha256":sha256_hex(&public_bytes),
        "checker_sha256":sha256_hex(&checker_bytes),
        "verified_controls":"separate field and curve arithmetic, all relation random streams, row ranks, scalar replay, source and destination transport, BSGS, rho returned scalars and caps"}));
}

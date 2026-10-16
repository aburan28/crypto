//! Independent arithmetic, inverse-fiber, and relation-witness replay.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::sync::OnceLock;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

type E = [u64; 2];
static CELL: OnceLock<(u64, u64, u64)> = OnceLock::new();
fn set_cell(p: u64, r: u64, h: u64) {
    assert!([
        (409, 433, 387),
        (421, 457, 387),
        (457, 421, 495),
        (1237, 1201, 1275),
        (2473, 2437, 2511),
        (5521, 5569, 5475),
        (9619, 9679, 9561),
        (19213, 19249, 19179)
    ]
    .contains(&(p, r, h)));
    assert!(CELL.set((p, r, h)).is_ok());
}
fn p() -> u64 {
    CELL.get().unwrap().0
}
fn r() -> u64 {
    CELL.get().unwrap().1
}
fn h() -> u64 {
    CELL.get().unwrap().2
}
const ARMS: [&str; 5] = [
    "inverse-velu",
    "full-image",
    "source-bsgs",
    "rho-negation",
    "rho-negation-aa",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Point {
    x: E,
    y: E,
    infinity: bool,
}

fn identity() -> Point {
    Point {
        x: [0, 0],
        y: [0, 0],
        infinity: true,
    }
}

fn parse(v: &Value) -> Point {
    if v["infinity"].as_bool().unwrap() {
        return identity();
    }
    Point {
        x: [v["x"][0].as_u64().unwrap(), v["x"][1].as_u64().unwrap()],
        y: [v["y"][0].as_u64().unwrap(), v["y"][1].as_u64().unwrap()],
        infinity: false,
    }
}

fn points(v: &Value) -> Vec<Point> {
    v.as_array().unwrap().iter().map(parse).collect()
}

fn tuple(q: Point) -> String {
    format!(
        "{},{},{},{},{}",
        u64::from(q.infinity),
        q.x[0],
        q.x[1],
        q.y[0],
        q.y[1]
    )
}

#[derive(Clone, Copy)]
struct Field {
    nu: u64,
}

impl Field {
    fn add(self, a: E, b: E) -> E {
        [(a[0] + b[0]) % p(), (a[1] + b[1]) % p()]
    }
    fn sub(self, a: E, b: E) -> E {
        [(a[0] + p() - b[0]) % p(), (a[1] + p() - b[1]) % p()]
    }
    fn mul(self, a: E, b: E) -> E {
        let p = p() as u128;
        [
            ((a[0] as u128 * b[0] as u128 + self.nu as u128 * a[1] as u128 * b[1] as u128) % p)
                as u64,
            ((a[0] as u128 * b[1] as u128 + a[1] as u128 * b[0] as u128) % p) as u64,
        ]
    }
    fn pow(self, mut a: u64, mut n: u64) -> u64 {
        let mut out = 1u64;
        while n > 0 {
            if n & 1 != 0 {
                out = out * a % p();
            }
            a = a * a % p();
            n >>= 1;
        }
        out
    }
    fn inv(self, a: E) -> E {
        let norm = (a[0] * a[0] + p() - self.nu * a[1] * a[1] % p()) % p();
        assert_ne!(norm, 0);
        self.mul([a[0], (p() - a[1]) % p()], [self.pow(norm, p() - 2), 0])
    }
    fn scale(self, a: E, k: u64) -> E {
        self.mul(a, [k % p(), 0])
    }
    fn neg(self, a: E) -> E {
        self.sub([0, 0], a)
    }
    fn pack(self, a: E) -> u64 {
        a[0] + p() * a[1]
    }
}

fn neg(f: Field, q: Point) -> Point {
    Point { y: f.neg(q.y), ..q }
}

fn add(f: Field, a: E, left: Point, right: Point) -> Point {
    if left.infinity {
        return right;
    }
    if right.infinity {
        return left;
    }
    if left.x == right.x && left.y == f.neg(right.y) {
        return identity();
    }
    let slope = if left == right {
        f.mul(
            f.add(f.scale(f.mul(left.x, left.x), 3), a),
            f.inv(f.scale(left.y, 2)),
        )
    } else {
        f.mul(f.sub(right.y, left.y), f.inv(f.sub(right.x, left.x)))
    };
    let x = f.sub(f.sub(f.mul(slope, slope), left.x), right.x);
    let y = f.sub(f.mul(slope, f.sub(left.x, x)), left.y);
    Point {
        x,
        y,
        infinity: false,
    }
}

fn scalar(f: Field, a: E, mut q: Point, mut k: u64) -> Point {
    let mut out = identity();
    while k > 0 {
        if k & 1 != 0 {
            out = add(f, a, out, q);
        }
        q = add(f, a, q, q);
        k >>= 1;
    }
    out
}

fn on_curve(f: Field, a: E, b: E, q: Point) -> bool {
    q.infinity || f.mul(q.y, q.y) == f.add(f.add(f.mul(f.mul(q.x, q.x), q.x), f.mul(a, q.x)), b)
}

fn velu(f: Field, a: E, t: Point, q: Point) -> Point {
    if q.infinity || q == t || q == neg(f, t) {
        return identity();
    }
    let plus = add(f, a, q, t);
    let minus = add(f, a, q, neg(f, t));
    assert!(!plus.infinity && !minus.infinity);
    Point {
        x: f.sub(f.add(f.add(q.x, plus.x), minus.x), f.scale(t.x, 2)),
        y: f.add(f.add(q.y, plus.y), minus.y),
        infinity: false,
    }
}

fn key(f: Field, q: Point) -> u64 {
    if q.infinity {
        return 0;
    }
    let sign = u64::from(f.pack(q.y) > f.pack(neg(f, q).y));
    ((f.pack(q.x) + 1) << 1) | sign
}

fn modmul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % r() as u128) as u64
}

fn pow_r(mut a: u64, mut n: u64) -> u64 {
    let mut out = 1u64;
    while n > 0 {
        if n & 1 != 0 {
            out = modmul(out, a);
        }
        a = modmul(a, a);
        n >>= 1;
    }
    out
}

fn solve_rows(rows: &[[u64; 10]]) -> Option<([u64; 9], usize)> {
    let mut matrix = rows.to_vec();
    let mut pivots = Vec::new();
    let mut lead = 0usize;
    for col in 0..9 {
        let Some(index) = (lead..matrix.len()).find(|&i| matrix[i][col] != 0) else {
            continue;
        };
        matrix.swap(lead, index);
        let inv = pow_r(matrix[lead][col], r() - 2);
        for j in col..10 {
            matrix[lead][j] = modmul(matrix[lead][j], inv);
        }
        for i in 0..matrix.len() {
            if i == lead || matrix[i][col] == 0 {
                continue;
            }
            let factor = matrix[i][col];
            for j in col..10 {
                matrix[i][j] = (matrix[i][j] + r() - modmul(factor, matrix[lead][j])) % r();
            }
        }
        pivots.push(col);
        lead += 1;
        if lead == 9 {
            break;
        }
    }
    assert!(
        matrix
            .iter()
            .skip(lead)
            .all(|row| row[..9].iter().any(|&v| v != 0) || row[9] == 0),
        "inconsistent relation"
    );
    if lead < 9 {
        return Some(([0; 9], lead));
    }
    let mut solution = [0u64; 9];
    for (i, &col) in pivots.iter().enumerate() {
        solution[col] = matrix[i][9];
    }
    Some((solution, lead))
}

fn base_for(public: &Value, record: &Value, arm: &str) -> Vec<Point> {
    if public["schema"] == "different_j_inverse_trend_public/v1" {
        let signed = points(&record["relation_setup"]["signed_points"]);
        assert_eq!(signed.len(), 16);
        let f = Field {
            nu: public["field_nu"].as_u64().unwrap(),
        };
        for i in 0..8 {
            assert_eq!(signed[2 * i + 1], neg(f, signed[2 * i]));
        }
        return signed.into_iter().step_by(2).collect();
    }
    let field = match arm {
        "source" | "isogeny-image" => "source_base_projected",
        "destination-direct" => "destination_base_projected",
        "full-image" | "inverse-velu" => "full_base_projected",
        _ => panic!("unknown relation arm"),
    };
    points(&public[field])
}

#[allow(clippy::too_many_arguments)]
fn check_relation(
    f: Field,
    public: &Value,
    record: &Value,
    row: &Value,
    arm: &str,
    g: Point,
    q: Point,
    mg: Point,
    mq: Point,
    t: Point,
) -> u64 {
    let source_a = [p() - 21, 0];
    let dest_a = [p() - 81, 0];
    let direct = ["destination-direct", "full-image", "inverse-velu"].contains(&arm);
    let query_a = if direct { dest_a } else { source_a };
    let table_a = if arm == "source" { source_a } else { dest_a };
    let query_g = if direct { mg } else { g };
    let query_q = if direct { mq } else { q };
    let base = base_for(public, record, arm);
    assert_eq!(base.len(), 8);
    let sign_a = if arm == "source" || arm == "isogeny-image" {
        source_a
    } else {
        dest_a
    };
    let signed: Vec<Point> = base.iter().flat_map(|&p| [p, neg(f, p)]).collect();
    for &point in &signed {
        assert!(scalar(f, sign_a, point, r()).infinity);
    }
    let table_points: Vec<Point> = if arm == "isogeny-image" {
        signed.iter().map(|&p| velu(f, source_a, t, p)).collect()
    } else {
        signed.clone()
    };
    let mut pair_entries = 0u64;
    let mut pair_keys = BTreeSet::new();
    let mut pair_table = BTreeMap::<u64, Vec<(usize, usize)>>::new();
    for i in 0..16 {
        for j in i..16 {
            let sum = add(f, table_a, table_points[i], table_points[j]);
            if !sum.infinity {
                pair_entries += 1;
                pair_keys.insert(key(f, sum));
                pair_table.entry(key(f, sum)).or_default().push((i, j));
            }
        }
    }
    assert_eq!(record["relation_setup"]["pair_entries"], pair_entries);
    assert_eq!(record["relation_setup"]["pair_keys"], pair_keys.len());
    assert_eq!(record["relation_setup"]["support_keys"], json!(pair_keys));
    let witnesses = row["relation_witnesses"].as_array().unwrap();
    assert_eq!(row["relations"], witnesses.len());
    assert_eq!(
        row["failed_queries"].as_u64().unwrap() + witnesses.len() as u64,
        row["trials"].as_u64().unwrap()
    );
    assert!(row["trials"].as_u64().unwrap() <= 8 * r());
    assert_eq!(row["invalid_candidates"], 0);
    let mut matrix = Vec::new();
    for witness in witnesses {
        let a = witness["a"].as_u64().unwrap();
        let b = witness["b"].as_u64().unwrap();
        let i = witness["i"].as_u64().unwrap() as usize;
        let j = witness["j"].as_u64().unwrap() as usize;
        assert!(a < r() && b < r() && i < 16 && j < 16 && i <= j);
        let query = add(
            f,
            query_a,
            scalar(f, query_a, query_g, a),
            scalar(f, query_a, query_q, b),
        );
        let projected = scalar(f, query_a, query, h());
        let searched = if arm == "isogeny-image" {
            velu(f, source_a, t, projected)
        } else {
            projected
        };
        let pair_sum = add(f, table_a, table_points[i], table_points[j]);
        assert_eq!(pair_sum, searched);
        if arm == "isogeny-image" {
            assert_eq!(add(f, source_a, signed[i], signed[j]), projected);
        }
        let mut lin = [0u64; 10];
        for index in [i, j] {
            lin[index / 2] = (lin[index / 2] + if index % 2 == 0 { 1 } else { r() - 1 }) % r();
        }
        lin[8] = (r() - modmul(h() % r(), b)) % r();
        lin[9] = modmul(h() % r(), a);
        matrix.push(lin);
    }
    let (solution, rank) = solve_rows(&matrix).unwrap();
    assert_eq!(rank, 9);
    assert_eq!(row["rank"], 9);
    assert_eq!(row["recovered"], solution[8]);
    if public["schema"] == "different_j_inverse_trend_public/v1" {
        let seed = row["seed"].as_u64().unwrap();
        let mut rng = StdRng::seed_from_u64(seed);
        let mut cursor = 0usize;
        for _ in 0..row["trials"].as_u64().unwrap() {
            let a = rng.gen_range(0..r());
            let b = rng.gen_range(0..r());
            let query = add(
                f,
                query_a,
                scalar(f, query_a, query_g, a),
                scalar(f, query_a, query_q, b),
            );
            let searched = scalar(f, query_a, query, h());
            if searched.infinity {
                continue;
            }
            if let Some(pairs) = pair_table.get(&key(f, searched)) {
                let &(i, j) = pairs.first().unwrap();
                let witness = &witnesses[cursor];
                assert_eq!(witness["a"], a);
                assert_eq!(witness["b"], b);
                assert_eq!(witness["i"], i);
                assert_eq!(witness["j"], j);
                assert_eq!(add(f, table_a, table_points[i], table_points[j]), searched);
                cursor += 1;
            }
        }
        assert_eq!(cursor, witnesses.len());
    }
    witnesses.len() as u64
}

fn check_bases(f: Field, public: &Value, t: Point) {
    for (raw_key, base_key, a, b) in [
        (
            "source_base_raw",
            "source_base_projected",
            [p() - 21, 0],
            [37, 0],
        ),
        (
            "destination_base_raw",
            "destination_base_projected",
            [p() - 81, 0],
            [p() - 243, 0],
        ),
        (
            "full_base_raw",
            "full_base_projected",
            [p() - 81, 0],
            [p() - 243, 0],
        ),
    ] {
        let raw = points(&public[raw_key]);
        let base = points(&public[base_key]);
        assert_eq!((raw.len(), base.len()), (8, 8));
        let mut orbits = BTreeSet::new();
        for (&p, &q) in raw.iter().zip(&base) {
            assert!(on_curve(f, a, b, p) && on_curve(f, a, b, q));
            assert_eq!(scalar(f, a, p, h()), q);
            assert!(!q.infinity && scalar(f, a, q, r()).infinity);
            assert!(orbits.insert(key(f, q).min(key(f, neg(f, q)))));
        }
    }
    let raw = points(&public["full_base_raw"]);
    let witnesses = points(&public["full_base_witness_source"]);
    assert_eq!(witnesses.len(), 8);
    for (&image, &witness) in raw.iter().zip(&witnesses) {
        assert!(on_curve(f, [p() - 21, 0], [37, 0], witness));
        assert_eq!(velu(f, [p() - 21, 0], t, witness), image);
        assert_ne!(witness.x, t.x);
        let x = witness.x;
        let u = image.x;
        let x2 = f.mul(x, x);
        let x3 = f.mul(x2, x);
        let cubic = f.sub(
            f.add(
                f.sub(x3, f.mul(f.add(u, [6, 0]), x2)),
                f.mul(f.add(f.scale(u, 6), [21, 0]), x),
            ),
            f.add(f.scale(u, 9), [32, 0]),
        );
        assert_eq!(cubic, [0, 0]);
        let mut orbit = [
            witness,
            add(f, [p() - 21, 0], witness, t),
            add(f, [p() - 21, 0], witness, neg(f, t)),
        ];
        assert!(orbit
            .iter()
            .all(|&point| !point.infinity && velu(f, [p() - 21, 0], t, point) == image));
        orbit.sort_by_key(|point| (point.x[0], point.x[1], point.y[0], point.y[1]));
        assert_eq!(orbit[0], witness);
    }
}

fn check_new_panel(public_path: &str, checker_path: &str, runs_path: &str) {
    let public_bytes = fs::read(public_path).unwrap();
    let checker_bytes = fs::read(checker_path).unwrap();
    let public: Value = serde_json::from_slice(&public_bytes).unwrap();
    let checker: Value = serde_json::from_slice(&checker_bytes).unwrap();
    assert_eq!(public["schema"], "different_j_inverse_trend_public/v1");
    assert_eq!(checker["schema"], "different_j_inverse_trend_checker/v1");
    assert_eq!(checker["public_sha256"], sha256_hex(&public_bytes));
    set_cell(
        public["p"].as_u64().unwrap(),
        public["r"].as_u64().unwrap(),
        public["cofactor"].as_u64().unwrap(),
    );
    let f = Field {
        nu: public["field_nu"].as_u64().unwrap(),
    };
    let source_a = [p() - 21, 0];
    let dest_a = [p() - 81, 0];
    let g = parse(&public["generator"]);
    let mg = parse(&public["mapped_generator"]);
    let t = parse(&public["kernel"][0]);
    assert!(on_curve(f, source_a, [37, 0], g));
    assert!(on_curve(f, dest_a, [p() - 243, 0], mg));
    assert_eq!(velu(f, source_a, t, g), mg);
    assert!(scalar(f, source_a, g, r()).infinity);
    assert!(scalar(f, dest_a, mg, r()).infinity);
    let targets: Vec<(Point, Point)> = public["targets"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| (parse(&row["source"]), parse(&row["destination"])))
        .collect();
    let answers: Vec<u64> = checker["scalars"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap())
        .collect();
    assert_eq!(targets.len(), 64);
    assert_eq!(answers.len(), 64);
    for (index, &(q, mq)) in targets.iter().enumerate() {
        assert_eq!(scalar(f, source_a, g, answers[index]), q);
        assert_eq!(scalar(f, dest_a, mg, answers[index]), mq);
        assert_eq!(velu(f, source_a, t, q), mq);
    }
    let root = std::path::Path::new(runs_path);
    let arms = [
        "inverse-velu",
        "source-bsgs",
        "rho-negation",
        "rho-negation-aa",
    ];
    let mut children = 0u64;
    let mut scalar_replays = 0u64;
    let mut relation_witnesses = 0u64;
    let mut censored = Vec::new();
    for round in 0u64..3 {
        for arm in arms {
            for prefix in [1, 64] {
                let id = format!("p{}_r{round}_{arm}_n{prefix}", p());
                let exit = fs::read_to_string(root.join(format!("{id}.exit"))).unwrap();
                if exit.trim() != "0" {
                    censored.push(json!({"id":id,"exit":exit.trim()}));
                    continue;
                }
                let record: Value =
                    serde_json::from_slice(&fs::read(root.join(format!("{id}.json"))).unwrap())
                        .unwrap();
                let meter: Value = serde_json::from_slice(
                    &fs::read(root.join(format!("{id}.meter.json"))).unwrap(),
                )
                .unwrap();
                assert_eq!(record["schema"], "different_j_trend_batch_run/v1");
                assert_eq!(record["p"], p());
                assert_eq!(record["r"], r());
                assert_eq!(record["h"], h());
                assert_eq!(record["public_sha256"], sha256_hex(&public_bytes));
                assert_eq!(record["round"], round);
                assert_eq!(record["arm"], arm);
                assert_eq!(record["prefix_length"], prefix);
                assert_eq!(record["relation_cap"], 8 * r());
                assert_eq!(record["rho_step_cap"], 20 * r());
                assert_eq!(record["all_verified"], true);
                assert_eq!(meter["exit_code"], 0);
                assert_eq!(meter["arm"], arm);
                assert_eq!(meter["round"], round);
                assert_eq!(meter["prefix_length"], prefix);
                assert!(meter["child_peak_rss_bytes"].as_u64().unwrap() <= 512 * 1024 * 1024);
                let rows = record["target_records"].as_array().unwrap();
                assert_eq!(rows.len(), prefix);
                if arm == "inverse-velu" {
                    let oracle = &record["relation_setup"]["oracle"];
                    let scanned = oracle["destination_points"].as_u64().unwrap();
                    let roots = oracle["rational_x_roots"].as_u64().unwrap();
                    let accepted = oracle["accepted_image_points"].as_u64().unwrap();
                    assert_eq!(oracle["polynomial_calls"], scanned);
                    assert_eq!(record["relation_setup"]["scan_count"], scanned);
                    assert_eq!(record["relation_setup"]["candidate_count"], accepted);
                    assert_eq!(
                        scanned,
                        accepted
                            + oracle["rejected_no_root"].as_u64().unwrap()
                            + oracle["rejected_no_lift"].as_u64().unwrap()
                            + oracle["rejected_map_mismatch"].as_u64().unwrap()
                    );
                    assert_eq!(
                        roots,
                        accepted
                            + oracle["rejected_no_lift"].as_u64().unwrap()
                            + oracle["rejected_map_mismatch"].as_u64().unwrap()
                    );
                }
                for (index, row) in rows.iter().enumerate() {
                    assert_eq!(row["position"], index);
                    let scalar_value = row["recovered"].as_u64().unwrap();
                    assert_eq!(scalar_value, answers[index]);
                    assert_eq!(scalar(f, source_a, g, scalar_value), targets[index].0);
                    assert_eq!(scalar(f, dest_a, mg, scalar_value), targets[index].1);
                    assert_eq!(row["verified"], true);
                    if arm == "inverse-velu" {
                        let seed = (0x2026_1015u64 ^ p())
                            .wrapping_add(round.wrapping_mul(0x9e37_79b9))
                            .wrapping_add((index as u64).wrapping_mul(0x0123_4567));
                        assert_eq!(row["seed"], seed);
                        relation_witnesses += check_relation(
                            f,
                            &public,
                            &record,
                            row,
                            arm,
                            g,
                            targets[index].0,
                            mg,
                            targets[index].1,
                            t,
                        );
                    } else if arm == "rho-negation" || arm == "rho-negation-aa" {
                        let attempts = row["attempts"].as_array().unwrap();
                        assert!(!attempts.is_empty() && attempts.len() <= 32);
                        assert_eq!(attempts.last().unwrap()["verified"], true);
                        for attempt in attempts {
                            assert!(attempt["steps"].as_u64().unwrap() <= 20 * r());
                        }
                    }
                    scalar_replays += 1;
                }
                children += 1;
            }
        }
    }
    println!(
        "{}",
        json!({"schema":"different_j_inverse_trend_run_check/v1",
        "status":if censored.is_empty(){"pass"}else{"partial"},
        "p":p(),"r":r(),"children":children,"scalar_replays":scalar_replays,
        "relation_witnesses":relation_witnesses,"censored":censored,
        "public_sha256":sha256_hex(&public_bytes),
        "checker_sha256":sha256_hex(&checker_bytes)})
    );
}

fn check_single_full_image(
    public_path: &str,
    checker_path: &str,
    run_path: &str,
    meter_path: &str,
) {
    let public_bytes = fs::read(public_path).unwrap();
    let checker_bytes = fs::read(checker_path).unwrap();
    let public: Value = serde_json::from_slice(&public_bytes).unwrap();
    let checker: Value = serde_json::from_slice(&checker_bytes).unwrap();
    let record: Value = serde_json::from_slice(&fs::read(run_path).unwrap()).unwrap();
    let meter: Value = serde_json::from_slice(&fs::read(meter_path).unwrap()).unwrap();
    assert_eq!(public["schema"], "different_j_inverse_trend_public/v1");
    assert_eq!(checker["schema"], "different_j_inverse_trend_checker/v1");
    assert_eq!(checker["public_sha256"], sha256_hex(&public_bytes));
    set_cell(
        public["p"].as_u64().unwrap(),
        public["r"].as_u64().unwrap(),
        public["cofactor"].as_u64().unwrap(),
    );
    let f = Field {
        nu: public["field_nu"].as_u64().unwrap(),
    };
    let source_a = [p() - 21, 0];
    let dest_a = [p() - 81, 0];
    let g = parse(&public["generator"]);
    let mg = parse(&public["mapped_generator"]);
    let q = parse(&public["targets"][0]["source"]);
    let mq = parse(&public["targets"][0]["destination"]);
    let t = parse(&public["kernel"][0]);
    assert!(on_curve(f, source_a, [37, 0], g));
    assert!(on_curve(f, dest_a, [p() - 243, 0], mg));
    assert_eq!(velu(f, source_a, t, g), mg);
    assert_eq!(velu(f, source_a, t, q), mq);
    assert_eq!(record["schema"], "different_j_trend_batch_run/v1");
    assert_eq!(record["public_sha256"], sha256_hex(&public_bytes));
    assert_eq!(record["arm"], "full-image");
    assert_eq!(record["round"], 0);
    assert_eq!(record["prefix_length"], 1);
    assert_eq!(record["all_verified"], true);
    assert_eq!(meter["schema"], "different_j_order_prefix_meter/v1");
    assert_eq!(meter["exit_code"], 0);
    assert_eq!(meter["arm"], "full-image");
    assert!(meter["child_peak_rss_bytes"].as_u64().unwrap() <= 512 * 1024 * 1024);
    let rows = record["target_records"].as_array().unwrap();
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row["position"], 0);
    assert_eq!(row["verified"], true);
    let scalar_value = row["recovered"].as_u64().unwrap();
    assert_eq!(scalar_value, checker["scalars"][0].as_u64().unwrap());
    assert_eq!(scalar(f, source_a, g, scalar_value), q);
    assert_eq!(scalar(f, dest_a, mg, scalar_value), mq);
    let witnesses = check_relation(f, &public, &record, row, "full-image", g, q, mg, mq, t);
    println!(
        "{}",
        json!({"schema":"different_j_full_image_single_check/v1",
        "status":"pass","p":p(),"r":r(),"scalar_replays":1,
        "relation_witnesses":witnesses,
        "public_sha256":sha256_hex(&public_bytes),
        "checker_sha256":sha256_hex(&checker_bytes)})
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("check-single-full-image") {
        assert_eq!(args.len(), 6);
        check_single_full_image(&args[2], &args[3], &args[4], &args[5]);
        return;
    }
    if args.get(1).map(String::as_str) == Some("check-new-panel") {
        assert_eq!(args.len(), 5);
        check_new_panel(&args[2], &args[3], &args[4]);
        return;
    }
    assert_eq!(
        args.len(),
        5,
        "usage: different_j_trend_batch_check PUBLIC CHECKER RUNS.jsonl METERS.jsonl"
    );
    let public_bytes = fs::read(&args[1]).unwrap();
    let checker_bytes = fs::read(&args[2]).unwrap();
    let run_bytes = fs::read(&args[3]).unwrap();
    let meter_bytes = fs::read(&args[4]).unwrap();
    let public: Value = serde_json::from_slice(&public_bytes).unwrap();
    let checker: Value = serde_json::from_slice(&checker_bytes).unwrap();
    assert_eq!(public["schema"], "different_j_trend_batch_public/v1");
    assert_eq!(checker["schema"], "different_j_trend_batch_checker/v1");
    set_cell(
        public["p"].as_u64().unwrap(),
        public["r"].as_u64().unwrap(),
        public["cofactor"].as_u64().unwrap(),
    );
    assert_eq!(checker["public_sha256"], sha256_hex(&public_bytes));
    assert_eq!(public["target_count"], 64);
    let nu = (2..p())
        .find(|&x| {
            let f = Field { nu: x };
            f.pow(x, (p() - 1) / 2) == p() - 1
        })
        .unwrap();
    let f = Field { nu };
    assert_eq!(public["field"], format!("F_{}[s]/(s^2-{nu})", p()));
    let g = parse(&public["generator"]);
    let mg = parse(&public["mapped_generator"]);
    let t = parse(&public["kernel"]);
    assert!(on_curve(f, [p() - 21, 0], [37, 0], g));
    assert!(on_curve(f, [p() - 81, 0], [p() - 243, 0], mg));
    assert!(
        scalar(f, [p() - 21, 0], g, r()).infinity && scalar(f, [p() - 81, 0], mg, r()).infinity
    );
    assert_eq!(velu(f, [p() - 21, 0], t, g), mg);
    check_bases(f, &public, t);
    let scalars: Vec<u64> = checker["scalars"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap())
        .collect();
    assert_eq!(scalars.len(), 64);
    let mut expected = Vec::new();
    for k in 1..r() {
        let q = scalar(f, [p() - 21, 0], g, k);
        let mq = scalar(f, [p() - 81, 0], mg, k);
        let key = tuple(q);
        let digest = sha256_hex(format!("different-j-prime-to-char-target-v2:13:{key}").as_bytes());
        expected.push((digest, key, q, mq, k));
    }
    expected.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
    let targets: Vec<(Point, Point)> = public["targets"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let q = parse(&row["source"]);
            let mq = parse(&row["destination"]);
            assert_eq!(
                (q, mq, scalars[i]),
                (expected[i].2, expected[i].3, expected[i].4)
            );
            assert_eq!(velu(f, [p() - 21, 0], t, q), mq);
            (q, mq)
        })
        .collect();
    let mut seen = BTreeSet::new();
    let mut scalar_checks = 0u64;
    let mut relation_checks = 0u64;
    let mut arm_counts = BTreeMap::<String, u64>::new();
    for line in run_bytes
        .split(|&c| c == b'\n')
        .filter(|line| !line.is_empty())
    {
        let record: Value = serde_json::from_slice(line).unwrap();
        assert_eq!(record["schema"], "different_j_trend_batch_run/v1");
        assert_eq!(record["construction"], "inverse-velu-v2");
        assert_eq!(
            (
                record["p"].as_u64(),
                record["r"].as_u64(),
                record["h"].as_u64()
            ),
            (Some(p()), Some(r()), Some(h()))
        );
        assert_eq!(record["relation_cap"], 8 * r());
        assert_eq!(record["rho_step_cap"], 20 * r());
        assert_eq!(record["public_sha256"], sha256_hex(&public_bytes));
        let round = record["round"].as_u64().unwrap();
        let arm = record["arm"].as_str().unwrap();
        let prefix = record["prefix_length"].as_u64().unwrap() as usize;
        assert!([1, 4, 16, 64].contains(&prefix));
        assert!(round < 3 && ARMS.contains(&arm));
        if arm == "inverse-velu" {
            let setup = &record["relation_setup"];
            let oracle = &setup["oracle"];
            let scanned = oracle["destination_points"].as_u64().unwrap();
            let roots = oracle["rational_x_roots"].as_u64().unwrap();
            let accepted = oracle["accepted_image_points"].as_u64().unwrap();
            assert_eq!(setup["scan_count"], scanned);
            assert_eq!(oracle["polynomial_calls"], scanned);
            assert_eq!(setup["candidate_count"], accepted);
            assert_eq!(accepted, public["full_candidate_count"].as_u64().unwrap());
            assert_eq!(
                scanned,
                accepted
                    + oracle["rejected_no_root"].as_u64().unwrap()
                    + oracle["rejected_no_lift"].as_u64().unwrap()
                    + oracle["rejected_map_mismatch"].as_u64().unwrap()
            );
            assert_eq!(
                roots,
                accepted
                    + oracle["rejected_no_lift"].as_u64().unwrap()
                    + oracle["rejected_map_mismatch"].as_u64().unwrap()
            );
            assert!(oracle["split_attempts"].as_u64().unwrap() <= 32 * scanned);
        }
        assert!(seen.insert((prefix, round, arm.to_string())));
        assert_eq!(record["target_count"], prefix);
        assert_eq!(record["all_verified"], true);
        let rows = record["target_records"].as_array().unwrap();
        assert_eq!(rows.len(), prefix);
        let prefixes = record["prefixes"].as_array().unwrap();
        let expected_prefixes: Vec<usize> = [1, 4, 16, 64]
            .into_iter()
            .filter(|&n| n <= prefix)
            .collect();
        assert_eq!(prefixes.len(), expected_prefixes.len());
        for (row, n) in prefixes.iter().zip(expected_prefixes) {
            assert_eq!(row["targets"], n);
            assert_eq!(row["verified"], true);
        }
        for (i, row) in rows.iter().enumerate() {
            assert_eq!(row["position"], i);
            let (q, mq) = targets[i];
            let k = row["recovered"].as_u64().unwrap();
            assert_eq!(k, scalars[i]);
            assert_eq!(scalar(f, [p() - 21, 0], g, k), q);
            assert_eq!(scalar(f, [p() - 81, 0], mg, k), mq);
            assert_eq!(row["verified"], true);
            if [
                "source",
                "isogeny-image",
                "destination-direct",
                "full-image",
                "inverse-velu",
            ]
            .contains(&arm)
            {
                let seed = (0x2026_1015u64 ^ p())
                    .wrapping_add(round.wrapping_mul(0x9e37_79b9))
                    .wrapping_add((i as u64).wrapping_mul(0x0123_4567));
                assert_eq!(row["seed"], seed);
                relation_checks += check_relation(f, &public, &record, row, arm, g, q, mg, mq, t);
            }
            scalar_checks += 1;
        }
        *arm_counts.entry(format!("{prefix}:{arm}")).or_default() += prefix as u64;
    }
    assert_eq!(seen.len(), 60);
    assert_eq!(scalar_checks, 1275);
    for prefix in [1, 4, 16, 64] {
        for arm in ARMS {
            assert_eq!(arm_counts[&format!("{prefix}:{arm}")], 3 * prefix as u64);
        }
    }
    let mut metered = BTreeSet::new();
    for line in meter_bytes
        .split(|&c| c == b'\n')
        .filter(|line| !line.is_empty())
    {
        let meter: Value = serde_json::from_slice(line).unwrap();
        assert_eq!(meter["schema"], "different_j_order_prefix_meter/v1");
        let prefix = meter["prefix_length"].as_u64().unwrap() as usize;
        let round = meter["round"].as_u64().unwrap();
        let arm = meter["arm"].as_str().unwrap().to_string();
        assert!(seen.contains(&(prefix, round, arm.clone())));
        assert!(metered.insert((prefix, round, arm)));
        assert_eq!(meter["exit_code"], 0);
        assert!(meter["child_cpu_us"].as_u64().unwrap() > 0);
        assert!(meter["stdout_bytes"].as_u64().unwrap() > 0);
        assert_eq!(meter["stderr_bytes"], 0);
        assert!(meter["child_peak_rss_bytes"].as_u64().unwrap() <= 512 * 1024 * 1024);
    }
    assert_eq!(metered, seen);
    println!(
        "{}",
        json!({"schema":"different_j_inverse_velu_check/v2","status":"pass","p":p(),"r":r(),"h":h(),
        "public_sha256":sha256_hex(&public_bytes),"checker_sha256":sha256_hex(&checker_bytes),
        "run_sha256":sha256_hex(&run_bytes),"meters_sha256":sha256_hex(&meter_bytes),
        "arm_rounds":seen.len(),"distinct_targets":targets.len(),
        "scalar_checks":scalar_checks,"relation_witness_checks":relation_checks,"arm_counts":arm_counts})
    );
}

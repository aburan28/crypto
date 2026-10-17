//! Standard reduced Tate/Miller diagnostic for the explicitly gated k=1,r=7
//! fixtures. Finite tests, not a new algebraic identity or general backend.
use crypto_lib::cryptanalysis::{
    curve_traits::arith,
    ecbench::canonical::sha256_hex,
    ic_boundary::{CountedGroup, GroupOps, IncrementalGauss, PrimeCurve, PrimePoint, RowStatus},
};
use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use std::{env, fs, path::Path, time::Instant};

#[derive(Default)]
struct Work {
    ec: GroupOps,
    muls: u64,
    squares: u64,
    inversions: u64,
    exponents: u64,
    adds_subs: u64,
    square_roots: u64,
    shift_candidates: u64,
    miller_steps: u64,
}
impl Work {
    fn mul(&mut self, a: u64, b: u64, p: u64) -> u64 {
        self.muls += 1;
        ((a as u128 * b as u128) % p as u128) as u64
    }
    fn sqr(&mut self, a: u64, p: u64) -> u64 {
        self.squares += 1;
        ((a as u128 * a as u128) % p as u128) as u64
    }
    fn add(&mut self, a: u64, b: u64, p: u64) -> u64 {
        self.adds_subs += 1;
        (a + b) % p
    }
    fn sub(&mut self, a: u64, b: u64, p: u64) -> u64 {
        self.adds_subs += 1;
        (a + p - b) % p
    }
    fn pow(&mut self, mut a: u64, mut n: u64, p: u64) -> u64 {
        self.exponents += 1;
        let mut b = 1;
        while n != 0 {
            if n & 1 != 0 {
                b = self.mul(b, a, p);
            }
            n >>= 1;
            if n != 0 {
                a = self.sqr(a, p);
            }
        }
        b
    }
    fn inv(&mut self, a: u64, p: u64) -> Option<u64> {
        if a == 0 {
            return None;
        }
        self.inversions += 1;
        Some(self.pow(a, p - 2, p))
    }
    fn point(&mut self, c: &PrimeCurve, x: u64) -> Option<PrimePoint> {
        let x2 = self.sqr(x, c.p);
        let x3 = self.mul(x2, x, c.p);
        let ax = self.mul(c.a, x, c.p);
        let t = self.add(x3, ax, c.p);
        let rhs = self.add(t, c.b, c.p);
        self.square_roots += 1;
        c.sqrt(rhs).map(|y| PrimePoint::affine(x, y))
    }
    fn json(&self) -> Value {
        json!({"ec_group":self.ec,"ec_field":self.ec.field_ops,
            "explicit_fp":{"multiplications":self.muls,"squarings":self.squares,
                "inversions":self.inversions,"exponentiations":self.exponents,
                "additions_subtractions":self.adds_subs},
            "explicit_fp_accounting":"Inversion exponentiation is expanded into the multiplication/square counts; do not price it twice. EC field ledgers are separate.",
            "sqrt_calls_unexpanded":self.square_roots,"shift_candidates":self.shift_candidates,
            "miller_steps":self.miller_steps})
    }
}

// g_(A,B)(R) = line_(A,B)(R) / vertical_(A+B)(R), with vertical special cases.
fn line(
    c: &PrimeCurve,
    a: PrimePoint,
    b: PrimePoint,
    q: PrimePoint,
    w: &mut Work,
) -> Option<(u64, PrimePoint)> {
    if q.infinity {
        return None;
    }
    if a.infinity {
        return Some((1, b));
    }
    if b.infinity {
        return Some((1, a));
    }
    let p = c.p;
    if a.x == b.x && (a.y != b.y || a.y == 0) {
        return Some((w.sub(q.x, a.x, p), PrimePoint::INFINITY));
    }
    let (num, den, sum) = if a == b {
        let x2 = w.sqr(a.x, p);
        let three_x2 = w.mul(3, x2, p);
        (
            w.add(three_x2, c.a, p),
            w.add(a.y, a.y, p),
            c.double(&mut w.ec, a),
        )
    } else {
        (
            w.sub(b.y, a.y, p),
            w.sub(b.x, a.x, p),
            c.add(&mut w.ec, a, b),
        )
    };
    let inv = w.inv(den, p)?;
    let slope = w.mul(num, inv, p);
    let dx = w.sub(q.x, a.x, p);
    let sx = w.mul(slope, dx, p);
    let dy = w.sub(q.y, a.y, p);
    let numerator = w.sub(dy, sx, p);
    if sum.infinity {
        return Some((numerator, sum));
    }
    let vertical = w.sub(q.x, sum.x, p);
    let inverse = w.inv(vertical, p)?;
    Some((w.mul(numerator, inverse, p), sum))
}

fn miller(c: &PrimeCurve, g: PrimePoint, q: PrimePoint, r: u64, w: &mut Work) -> Option<u64> {
    let mut f = 1;
    let mut t = g;
    let top = 63 - r.leading_zeros();
    for bit in (0..top).rev() {
        w.miller_steps += 1;
        let (v, next) = line(c, t, t, q, w)?;
        let square = w.sqr(f, c.p);
        f = w.mul(square, v, c.p);
        t = next;
        if r & (1 << bit) != 0 {
            let (v, next) = line(c, t, g, q, w)?;
            f = w.mul(f, v, c.p);
            t = next;
        }
    }
    if !t.infinity || f == 0 {
        None
    } else {
        Some(f)
    }
}

fn tate(c: &PrimeCurve, g: PrimePoint, q: PrimePoint, r: u64, w: &mut Work) -> Result<u64, String> {
    if (c.p - 1) % r != 0 || q.infinity {
        return Err("k=1 requires r | p-1 and nonzero target".into());
    }
    for x in 1..=128 {
        w.shift_candidates += 1;
        let Some(s) = w.point(c, x) else {
            continue;
        };
        let qs = c.add(&mut w.ec, q, s);
        let Some(num) = miller(c, g, qs, r, w) else {
            continue;
        };
        let Some(den) = miller(c, g, s, r, w) else {
            continue;
        };
        let Some(inv) = w.inv(den, c.p) else {
            continue;
        };
        let ratio = w.mul(num, inv, c.p);
        let value = w.pow(ratio, (c.p - 1) / r, c.p);
        if value != 0 {
            return Ok(value);
        }
    }
    Err("shift search exhausted".into())
}

fn projected_point(c: &PrimeCurve, h: u64, start: u64, w: &mut Work) -> Result<PrimePoint, String> {
    for x in start..start + 10000 {
        if let Some(p) = w.point(c, x) {
            let q = c.mul(&mut w.ec, p, h);
            if !q.infinity {
                return Ok(q);
            }
        }
    }
    Err("point projection budget exhausted".into())
}

fn relation_solve(
    p: u64,
    r: u64,
    alpha: u64,
    beta: u64,
    w: &mut Work,
) -> Result<(u64, Value), String> {
    if r != 7 || alpha == 1 || w.pow(alpha, r, p) != 1 || w.pow(beta, r, p) != 1 {
        return Err("diagnostic requires nontrivial generator in order-7 destination".into());
    }
    let mut basis = vec![];
    let mut scanned = 0;
    for prime in [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97,
    ] {
        scanned += 1;
        let z = w.pow(prime, (p - 1) / r, p);
        if z == 1 {
            continue;
        }
        let inverse = w.inv(z, p).unwrap();
        if basis
            .iter()
            .any(|(_, v): &(u64, u64)| *v == z || *v == inverse)
        {
            continue;
        }
        basis.push((prime, z));
        if basis.len() == 2 {
            break;
        }
    }
    if basis.len() != 2 {
        return Err("factor-base scan exhausted".into());
    }
    let mut signed = vec![];
    for (i, (_, z)) in basis.iter().enumerate() {
        signed.push((*z, i, 1));
        signed.push((w.inv(*z, p).unwrap(), i, r - 1));
    }
    let mut gauss = IncrementalGauss::new(3, r);
    let mut relations = vec![];
    let mut calls = 0;
    let mut pair_tests = 0;
    let mut hits = 0;
    let mut identity_skips = 0;
    let mut walk = beta;
    'queries: for shift in 0..r {
        if walk == 1 {
            identity_skips += 1;
        } else {
            calls += 1;
            let mut hit = false;
            for &(left, i, e) in &signed {
                for &(right, j, f) in &signed {
                    pair_tests += 1;
                    if w.mul(left, right, p) != walk {
                        continue;
                    }
                    hit = true;
                    let mut row = vec![0; 3];
                    row[i] = (row[i] + e) % r;
                    row[j] = (row[j] + f) % r;
                    row[2] = r - 1;
                    relations.push(
                        json!({"shift":shift,"row":row,"rhs":shift,"destination_value":walk}),
                    );
                    if gauss.add_row(row, shift) == RowStatus::Inconsistent {
                        return Err("inconsistent checked relations".into());
                    }
                    if gauss.rank() == 3 {
                        hits += 1;
                        break 'queries;
                    }
                }
            }
            if hit {
                hits += 1;
            }
        }
        walk = w.mul(walk, alpha, p);
    }
    if gauss.rank() != 3 {
        return Err(format!("rank exhausted: {}", gauss.rank()));
    }
    let x = gauss.pinned(2).ok_or("target column unpinned")?;
    for (i, (_, z)) in basis.iter().enumerate() {
        let log = gauss.pinned(i).ok_or("factor-base column unpinned")?;
        if w.pow(alpha, log, p) != *z {
            return Err("factor-base log failed verification".into());
        }
    }
    if w.pow(alpha, x, p) != beta {
        return Err("destination answer failed verification".into());
    }
    Ok((
        x,
        json!({"factor_base":basis,"factor_base_size":2,"prime_candidates_scanned":scanned,
        "queries":calls,"queries_with_relation":hits,"candidate_pair_tests":pair_tests,
            "relation_count":relations.len(),"relations_per_query":relations.len() as f64 / calls as f64,
            "query_hit_yield":hits as f64 / calls as f64,
            "relations":relations,"rank":gauss.rank(),"columns":3,"matrix_entry_ops":gauss.row_ops,
        "dependent_rows":gauss.dependent,"identity_queries_skipped":identity_skips,
        "factor_base_logs_checked":true,"destination_answer_checked":true,
        "yield_definition":"relations.len / queries; multiple relations per query allowed; also report queries_with_relation / queries"}),
    ))
}

fn run(c: &PrimeCurve, n: u64, target_index: u64) -> Result<Value, String> {
    let start = Instant::now();
    let factors = arith::factor(&n.into(), 1_000_000);
    if !factors.complete()
        || factors.largest_prime().unwrap().to_string() != "7"
        || n % 7 != 0
        || (n / 7) % 7 == 0
    {
        return Err("only certified largest-prime order-7 fixtures supported".into());
    }
    let mut setup = Work::default();
    let g = projected_point(c, n / 7, 1, &mut setup)?;
    let q = projected_point(c, n / 7, 10001 + 10000 * target_index, &mut setup)?;
    if !c.mul(&mut setup.ec, g, 7).infinity || !c.mul(&mut setup.ec, q, 7).infinity {
        return Err("subgroup certificate failed".into());
    }
    let setup_ns = start.elapsed().as_nanos().to_string();
    let mut pairing = Work::default();
    let pair_start = Instant::now();
    let alpha = tate(c, g, g, 7, &mut pairing)?;
    let beta = tate(c, g, q, 7, &mut pairing)?;
    let pair_ns = pair_start.elapsed().as_nanos().to_string();
    let mut ff = Work::default();
    let ff_start = Instant::now();
    let (x, relations) = relation_solve(c.p, 7, alpha, beta, &mut ff)?;
    let ff_ns = ff_start.elapsed().as_nanos().to_string();
    let mut verification = Work::default();
    let verify_start = Instant::now();
    if c.mul(&mut verification.ec, g, x) != q {
        return Err("recovered scalar failed EC verification".into());
    }
    let verify_ns = verify_start.elapsed().as_nanos().to_string();
    Ok(
        json!({"curve":c,"order":n,"r":7,"generator":[g.x,g.y],"target":[q.x,q.y],
        "target_kind":"public deterministic point search and cofactor projection; no scalar supplied to solver",
        "target_index":target_index,"alpha":alpha,"beta":beta,"recovered":x,"status":"verified",
        "phases":{"setup":{"operations":setup.json(),"elapsed_ns":setup_ns},
            "pairing":{"operations":pairing.json(),"elapsed_ns":pair_ns},
            "factor_base_relations_matrix_and_ff_verification":{"operations":ff.json(),"elapsed_ns":ff_ns},
            "ec_verification":{"operations":verification.json(),"elapsed_ns":verify_ns}},
        "end_to_end_ns":start.elapsed().as_nanos().to_string(),"relation_report":relations,
        "unpriced":["native input order factorization", "square-root internals (calls counted)", "allocation and integer comparisons", "matrix inversions not expanded; matrix entry operations counted"],
        "timing_status":"UNQUALIFIED_DESCRIPTIVE; not an ecbench promotion or external speed comparison"}),
    )
}

fn main() {
    let args: Vec<_> = env::args().collect();
    if args.get(1).map(String::as_str) == Some("compare") {
        assert_eq!(
            args.len(),
            5,
            "finite-field-k1 compare BASE CANDIDATE FRESH_RECEIPT"
        );
        assert!(!Path::new(&args[4]).exists());
        let raw_a = fs::read(&args[2]).unwrap();
        let raw_b = fs::read(&args[3]).unwrap();
        let mut a: Value = serde_json::from_slice(&raw_a).unwrap();
        let mut b: Value = serde_json::from_slice(&raw_b).unwrap();
        let records = a["rows"].as_array().unwrap().len();
        fn omit_times(v: &mut Value) {
            match v {
                Value::Object(m) => {
                    m.remove("end_to_end_ns");
                    m.remove("elapsed_ns");
                    for v in m.values_mut() {
                        omit_times(v);
                    }
                }
                Value::Array(a) => {
                    for v in a {
                        omit_times(v);
                    }
                }
                _ => {}
            }
        }
        omit_times(&mut a);
        omit_times(&mut b);
        let ok = a == b;
        let receipt = json!({"schema":"volcano-ff-k1-replay-match/v1","ok":ok,
            "base_sha256":sha256_hex(&raw_a),"candidate_sha256":sha256_hex(&raw_b),"records":records,
            "scope":"Exact JSON equality after excluding elapsed_ns and end_to_end_ns at every level; includes identities, outcomes, relations and all operation counters. Same binary, not independent implementation."});
        fs::write(&args[4], serde_json::to_vec_pretty(&receipt).unwrap()).unwrap();
        assert!(ok, "re-execution differs; failed receipt preserved");
        return;
    }
    assert!(
        args.len() == 3 || args.len() == 4,
        "finite-field-k1 ORIGINAL_C_JSON FRESH_OUTPUT [FROZEN_C_PREFIX_DB]"
    );
    let out = Path::new(&args[2]);
    assert!(!out.exists());
    let raw = fs::read(&args[1]).unwrap();
    let v: Value = serde_json::from_slice(&raw).unwrap();
    let c = PrimeCurve {
        p: v["p"].as_u64().unwrap(),
        a: v["a"].as_u64().unwrap(),
        b: v["b"].as_u64().unwrap(),
    };
    assert_eq!(
        (c.p, c.a, c.b),
        (1637958981148213249, 571607911220311176, 29329524900995129)
    );
    let n = v["order"].as_u64().unwrap();
    assert_eq!(n, 1637958981148213248);
    let mut cases = vec![(None, c, false, None)];
    let prefix = if args.len() == 4 {
        let raw_db = fs::read(&args[3]).unwrap();
        let digest = sha256_hex(&raw_db);
        assert_eq!(
            digest,
            "468ce8320b6b7deae37198ca3afc07534fa4fe0e950f621e45ad322f23a2c881"
        );
        let uri = format!(
            "file:{}?immutable=1",
            Path::new(&args[3]).canonicalize().unwrap().display()
        );
        let db = Connection::open_with_flags(
            uri,
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_URI,
        )
        .unwrap();
        let meta = |key: &str| -> Value {
            let s: String = db
                .query_row("SELECT value FROM metadata WHERE key=?", [key], |r| {
                    r.get(0)
                })
                .unwrap();
            serde_json::from_str(&s).unwrap()
        };
        let plan = meta("plan");
        let receipt = meta("verification");
        assert_eq!(plan["input"], v);
        assert_eq!(plan["order"], n);
        assert_eq!(receipt["status"], "REPLAYED");
        assert_eq!(receipt["vertex_class_complete"], false);
        let mut stmt = db
            .prepare("SELECT id,a,b,expanded,conductor FROM nodes ORDER BY id")
            .unwrap();
        cases = stmt
            .query_map([], |r| {
                Ok((
                    Some(r.get::<_, u64>(0)?),
                    PrimeCurve {
                        p: 1637958981148213249,
                        a: r.get(1)?,
                        b: r.get(2)?,
                    },
                    r.get::<_, u64>(3)? == 1,
                    r.get::<_, Option<u64>>(4)?,
                ))
            })
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(cases.len(), 1869);
        Some(
            json!({"sha256":digest,"verification":receipt,"vertices":cases.len(),
            "scope":"Every stored prefix vertex; original C frontier is NOT exhausted; conductor strata NOT_ESTABLISHED"}),
        )
    } else {
        None
    };
    let requested = cases.len() * 2;
    let mut rows = vec![];
    for (vertex, c, expanded, conductor) in cases {
        for i in 0..2 {
            let mut result = match run(&c, n, i) {
                Ok(v) => v,
                Err(e) => json!({"curve":c,"target_index":i,"status":"failed","reason":e}),
            };
            result["graph_vertex"] = json!(vertex);
            result["prefix_expanded"] = json!(expanded);
            result["raw_prefix_conductor_uncertified"] = json!(conductor);
            rows.push(result);
        }
    }
    // Bind the source compiled into this executable, not a mutable runtime file.
    let source = include_bytes!("finite_field_k1.rs");
    let verified: Vec<_> = rows.iter().filter(|v| v["status"] == "verified").collect();
    let stats = |values: Vec<u64>| -> Value {
        if values.is_empty() {
            return json!({"n":0,"min":null,"max":null,"mean":null});
        }
        let sum: u128 = values.iter().map(|v| *v as u128).sum();
        json!({"n":values.len(),"min":values.iter().min(),"max":values.iter().max(),
            "sum":sum.to_string(),"mean":sum as f64 / values.len() as f64})
    };
    let mut metrics = serde_json::Map::new();
    for key in [
        "relation_count",
        "queries",
        "candidate_pair_tests",
        "matrix_entry_ops",
    ] {
        metrics.insert(
            key.into(),
            stats(
                verified
                    .iter()
                    .map(|v| v["relation_report"][key].as_u64().unwrap())
                    .collect(),
            ),
        );
    }
    for key in ["relations_per_query", "query_hit_yield"] {
        let values: Vec<f64> = verified
            .iter()
            .map(|v| v["relation_report"][key].as_f64().unwrap())
            .collect();
        let summary = if values.is_empty() {
            json!({"n":0,"min":null,"max":null,"mean":null})
        } else {
            json!({"n":values.len(),"min":values.iter().copied().reduce(f64::min),
                "max":values.iter().copied().reduce(f64::max),"mean":values.iter().sum::<f64>()/values.len() as f64})
        };
        metrics.insert(key.into(), summary);
    }
    let mut phase_statistics = serde_json::Map::new();
    for phase in [
        "setup",
        "pairing",
        "factor_base_relations_matrix_and_ff_verification",
        "ec_verification",
    ] {
        for field in [
            "ec_group/adds",
            "ec_group/doubles",
            "ec_group/scalar_mults",
            "ec_field/muls",
            "ec_field/sqrs",
            "ec_field/invs",
            "explicit_fp/multiplications",
            "explicit_fp/squarings",
            "explicit_fp/inversions",
            "explicit_fp/exponentiations",
            "explicit_fp/additions_subtractions",
            "sqrt_calls_unexpanded",
            "shift_candidates",
            "miller_steps",
        ] {
            let pointer = format!("/phases/{phase}/operations/{field}");
            phase_statistics.insert(
                format!("{phase}/{field}"),
                stats(
                    verified
                        .iter()
                        .map(|v| v.pointer(&pointer).unwrap().as_u64().unwrap())
                        .collect(),
                ),
            );
        }
    }
    let times = stats(
        verified
            .iter()
            .map(|v| v["end_to_end_ns"].as_str().unwrap().parse().unwrap())
            .collect(),
    );
    fs::write(out,serde_json::to_vec_pretty(&json!({"schema":"volcano-ff-k1-diagnostic/v1",
        "input_sha256":sha256_hex(&raw),"source_sha256":sha256_hex(source),
        "scope":if prefix.is_some(){"All 1869 known original C prefix vertices; not its full graph or all conductor strata"}else{"Original C source model only; not its full graph or all conductor strata"},
        "prefix":prefix,
        "backend":"Reduced Tate/Miller then projected-prime-subgroup two-summand finite-field relations; no generic solver fallback",
        "summary":{"requested":requested,"verified":verified.len(),"failed":requested-verified.len(),
            "end_to_end_ns":times,"metrics":metrics,"phase_operations":phase_statistics,"scope":"Present verified observations only; raw failed rows preserved. Two public target searches per model; not full-class or qualified timing."},
        "rows":rows})).unwrap()).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    fn all_scalars(c: PrimeCurve, n: u64) {
        let mut setup = Work::default();
        let g = projected_point(&c, n / 7, 1, &mut setup).unwrap();
        let mut w = Work::default();
        let alpha = tate(&c, g, g, 7, &mut w).unwrap();
        assert_ne!(alpha, 1);
        for x in 1..7 {
            let q = c.mul(&mut w.ec, g, x);
            let beta = tate(&c, g, q, 7, &mut w).unwrap();
            assert_eq!(
                w.pow(alpha, x, c.p),
                beta,
                "finite fixture pairing evaluation"
            );
            let (found, receipt) = relation_solve(c.p, 7, alpha, beta, &mut w).unwrap();
            assert_eq!(found, x);
            assert_eq!(receipt["rank"], 3);
            assert_eq!(c.mul(&mut w.ec, g, found), q);
        }
    }
    #[test]
    fn control_all_nonzero_inputs() {
        all_scalars(
            PrimeCurve {
                p: 113,
                a: 19,
                b: 2,
            },
            112,
        );
    }
    #[test]
    fn original_c_all_nonzero_inputs() {
        all_scalars(
            PrimeCurve {
                p: 1637958981148213249,
                a: 571607911220311176,
                b: 29329524900995129,
            },
            1637958981148213248,
        );
    }
    #[test]
    fn rejects_degenerate_destination() {
        assert!(relation_solve(113, 7, 1, 1, &mut Work::default()).is_err());
    }
    #[test]
    fn destination_relations_check_every_nontrivial_generator_and_target() {
        for (c, n) in [
            (
                PrimeCurve {
                    p: 113,
                    a: 19,
                    b: 2,
                },
                112,
            ),
            (
                PrimeCurve {
                    p: 1637958981148213249,
                    a: 571607911220311176,
                    b: 29329524900995129,
                },
                1637958981148213248,
            ),
        ] {
            let mut w = Work::default();
            let g = projected_point(&c, n / 7, 1, &mut w).unwrap();
            let base = tate(&c, g, g, 7, &mut w).unwrap();
            for generator_power in 1..7 {
                let alpha = w.pow(base, generator_power, c.p);
                for x in 1..7 {
                    let beta = w.pow(alpha, x, c.p);
                    assert_eq!(relation_solve(c.p, 7, alpha, beta, &mut w).unwrap().0, x);
                }
            }
        }
    }
}

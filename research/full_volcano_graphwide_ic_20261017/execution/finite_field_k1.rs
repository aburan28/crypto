//! Standard reduced Tate/Miller diagnostic for the explicitly gated k=1,r=7
//! fixtures. Finite tests, not a new algebraic identity or general backend.
use crypto_lib::cryptanalysis::{
    curve_traits::arith,
    ecbench::canonical::sha256_hex,
    ic_boundary::{CountedGroup, GroupOps, IncrementalGauss, PrimeCurve, PrimePoint, RowStatus},
};
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
    assert_eq!(
        args.len(),
        3,
        "finite-field-k1 ORIGINAL_C_JSON FRESH_OUTPUT"
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
    let mut rows = vec![];
    for i in 0..2 {
        rows.push(match run(&c, n, i) {
            Ok(v) => v,
            Err(e) => json!({"target_index":i,"status":"failed","reason":e}),
        });
    }
    let source = fs::read(file!()).expect("run from source checkout so source is bound");
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
    let times = stats(
        verified
            .iter()
            .map(|v| v["end_to_end_ns"].as_str().unwrap().parse().unwrap())
            .collect(),
    );
    fs::write(out,serde_json::to_vec_pretty(&json!({"schema":"volcano-ff-k1-diagnostic/v1",
        "input_sha256":sha256_hex(&raw),"source_sha256":sha256_hex(&source),
        "scope":"Original C source model only; not its full graph or all conductor strata",
        "backend":"Reduced Tate/Miller then projected-prime-subgroup two-summand finite-field relations; no generic solver fallback",
        "summary":{"requested":2,"verified":verified.len(),"failed":2-verified.len(),
            "end_to_end_ns":times,"metrics":metrics,"scope":"Present verified observations only; raw failed rows preserved. Two public targets on one original model, not graph-wide or qualified timing."},
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
}

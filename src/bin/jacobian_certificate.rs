//! Exact, bounded Hermitian-lattice certificates; geometric conclusions stay conditional.
use clap::Parser;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{collections::HashSet, fs, path::PathBuf};

type Check<T> = Result<T, String>;
type Element = [i128; 2];

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Binding {
    slug: String,
    curve_uid: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Candidate {
    label: String,
    discriminant: i64,
    a: i64,
    d: i64,
    /// b = b[0] + b[1]*omega; omega=(t+sqrt(D))/2, t=D mod 2.
    b: [i64; 2],
    binding: Option<Binding>,
    /// References are retained as unverified provenance, never as proof flags.
    evidence_refs: Vec<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    schema_version: String,
    candidates: Vec<Candidate>,
}

fn digest(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

fn sqrt(n: i128) -> i128 {
    let (mut lo, mut hi) = (0, n + 1);
    while hi - lo > 1 {
        let m = (lo + hi) / 2;
        if m <= n / m {
            lo = m;
        } else {
            hi = m;
        }
    }
    lo
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

fn squarefree(n: i128) -> bool {
    let mut p = 2;
    while p * p <= n {
        if n % (p * p) == 0 {
            return false;
        }
        p += 1;
    }
    true
}

fn fundamental(d: i128) -> bool {
    if d >= 0 {
        return false;
    }
    if d.rem_euclid(4) == 1 {
        squarefree(-d)
    } else if d % 4 == 0 {
        let m = d / 4;
        matches!(m.rem_euclid(4), 2 | 3) && squarefree(-m)
    } else {
        false
    }
}

fn reduced_forms(d: i128) -> Vec<[i128; 3]> {
    let mut forms = vec![];
    for a in 1..=sqrt((-d) / 3) {
        for b in -a..=a {
            if (b * b - d) % (4 * a) != 0 {
                continue;
            }
            let c = (b * b - d) / (4 * a);
            if c < a || gcd(gcd(a, b), c) != 1 {
                continue;
            }
            if (b.abs() == a || a == c) && b < 0 {
                continue;
            }
            forms.push([a, b, c]);
        }
    }
    forms
}

fn norm(x: Element, t: i128, n: i128) -> i128 {
    x[0] * x[0] + t * x[0] * x[1] + n * x[1] * x[1]
}

fn mul(x: Element, y: Element, t: i128, n: i128) -> Element {
    [
        x[0] * y[0] - n * x[1] * y[1],
        x[0] * y[1] + x[1] * y[0] + t * x[1] * y[1],
    ]
}

fn charge(work: &mut u64, budget: u64) -> Check<()> {
    if *work >= budget {
        return Err("exhaustive enumeration budget exhausted".into());
    }
    *work += 1;
    Ok(())
}

fn elements(limit: i128, d: i128, work: &mut u64, budget: u64) -> Check<Vec<(Element, i128)>> {
    let t = d.rem_euclid(2);
    let n = (t * t - d) / 4;
    let vb = sqrt(4 * limit / -d);
    let mut out = vec![];
    // 4*N(u+v*w)=(2u+t*v)^2+|D|*v^2 bounds both coordinates exactly.
    for v in -vb..=vb {
        let ub = sqrt(limit) + (v.abs() + 1) / 2 + 1;
        for u in -ub..=ub {
            charge(work, budget)?;
            let q = norm([u, v], t, n);
            if q <= limit {
                out.push(([u, v], q));
            }
        }
    }
    Ok(out)
}

fn minimum(c: &Candidate, work: &mut u64, budget: u64) -> Check<(i128, [Element; 2])> {
    let (a, d) = (i128::from(c.a), i128::from(c.discriminant));
    let t = d.rem_euclid(2);
    let n = (t * t - d) / 4;
    let b = c.b.map(i128::from);
    let mut best = i128::from(c.a.min(c.d));
    let mut witness = if c.a <= c.d {
        [[1, 0], [0, 0]]
    } else {
        [[0, 0], [1, 0]]
    };
    let candidates = elements(a * best, d, work, budget)?;
    // det(H)=1 gives a*h(x,y)=N(z)+N(y), z=a*x+b*y.
    // Enumerating both norm balls and the congruence z=b*y mod a is complete.
    for &(y, ny) in &candidates {
        let by = mul(b, y, t, n);
        for &(z, nz) in &candidates {
            charge(work, budget)?;
            let total = ny + nz;
            if total == 0 || total >= a * best {
                continue;
            }
            let numerator = [z[0] - by[0], z[1] - by[1]];
            if numerator.iter().any(|q| q % a != 0) {
                continue;
            }
            if total % a != 0 {
                return Err("internal norm identity failure".into());
            }
            best = total / a;
            witness = [[numerator[0] / a, numerator[1] / a], y];
        }
    }
    Ok((best, witness))
}

fn analyze(c: &Candidate, budget: u64) -> Value {
    let mut result = json!({"status":"unsupported", "principal":null,
        "minimum":null,"minimum_witness":null,"geometrically_indecomposable":null,
        "reduced_forms":null,"decomposable_minimum_upper_bound":null,"work_units":0});
    let d = i128::from(c.discriminant);
    // These public input bounds keep every intermediate comfortably inside i128.
    if !(-1_000_000..=-3).contains(&d)
        || !(1..=1000).contains(&c.a)
        || !(1..=1000).contains(&c.d)
        || c.b.iter().any(|b| !(-1000..=1000).contains(b))
    {
        result["reason"] =
            json!("supported bounds: -1000000 <= D <= -3; 1 <= a,d <= 1000; |b_i| <= 1000");
        return result;
    }
    if !fundamental(d) {
        result["reason"] = json!("requires a negative fundamental discriminant (maximal order)");
        return result;
    }
    let t = d.rem_euclid(2);
    let n = (t * t - d) / 4;
    let determinant = i128::from(c.a) * i128::from(c.d) - norm(c.b.map(i128::from), t, n);
    result["determinant"] = json!(determinant);
    result["principal"] = json!(determinant == 1);
    if determinant != 1 {
        result["status"] = json!("not_principal_candidate");
        result["reason"] = json!("candidate does not have determinant one");
        return result;
    }
    let forms = reduced_forms(d);
    let bound = forms
        .iter()
        .map(|f| f[0])
        .max()
        .expect("principal ideal class exists");
    result["reduced_forms"] = json!(forms);
    result["decomposable_minimum_upper_bound"] = json!(bound);
    let mut work = 0;
    match minimum(c, &mut work, budget) {
        Ok((min, witness)) => {
            result["minimum"] = json!(min);
            result["minimum_witness"] = json!(witness);
            let (status, indecomposable) = if min > bound {
                ("indecomposable_lattice", Some(true))
            } else if min == 1 {
                ("decomposable_lattice", Some(false))
            } else {
                ("inconclusive", None)
            };
            result["status"] = json!(status);
            result["geometrically_indecomposable"] = json!(indecomposable);
        }
        Err(reason) => {
            result["status"] = json!("inconclusive_budget");
            result["reason"] = json!(reason);
        }
    }
    result["work_units"] = json!(work);
    result
}

fn binding(c: &Candidate, registry: Option<&Value>) -> Check<Value> {
    let Some(b) = &c.binding else {
        return Ok(Value::Null);
    };
    let registry = registry.ok_or("a binding requires --registry")?;
    if registry["schema_version"] != 1 {
        return Err("unsupported registry version".into());
    }
    let rows = registry["curves"]
        .as_array()
        .ok_or("missing registry curves")?;
    let rows: Vec<_> = rows.iter().filter(|r| r["slug"] == b.slug).collect();
    if rows.len() != 1 {
        return Err("curve slug missing or ambiguous".into());
    }
    let row = rows[0];
    let reps = row["representations"]
        .as_array()
        .ok_or("missing representations")?;
    let reps: Vec<_> = reps
        .iter()
        .filter(|r| r["curve_uid"] == b.curve_uid)
        .collect();
    if reps.len() != 1 {
        return Err("full curve UID missing or ambiguous for slug".into());
    }
    let rep = reps[0];
    let uid = b
        .curve_uid
        .strip_prefix("urn:ec-record:1:sha256:")
        .ok_or("invalid curve UID")?;
    if uid.len() != 64 || !uid.bytes().all(|x| x.is_ascii_hexdigit()) {
        return Err("invalid curve UID digest".into());
    }
    if !rep["ec1"].is_string() || !row["model_json"].is_string() {
        return Err("missing exact catalog identity metadata".into());
    }
    Ok(
        json!({"slug":b.slug,"curve_uid":b.curve_uid,"ec1":rep["ec1"],
        "model_json":row["model_json"],"field":rep["field"],"curve":rep["curve"],
        "status":"matched_catalog_metadata_only"}),
    )
}

fn report(input: &[u8], registry: Option<&[u8]>, budget: u64) -> Check<Value> {
    if !(1..=50_000_000).contains(&budget) {
        return Err("budget must be 1..50000000".into());
    }
    let parsed: Input = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if parsed.schema_version != "jacobian-candidates/v1"
        || parsed.candidates.is_empty()
        || parsed.candidates.len() > 100
    {
        return Err("expected jacobian-candidates/v1 with 1..100 candidates".into());
    }
    let reg: Option<Value> = registry
        .map(serde_json::from_slice)
        .transpose()
        .map_err(|e| e.to_string())?;
    let mut seen = HashSet::new();
    let mut records = vec![];
    for c in parsed.candidates {
        if c.label.is_empty() || !seen.insert(c.label.clone()) {
            return Err("candidate labels must be nonempty and unique".into());
        }
        let lattice = analyze(&c, budget);
        let status = if lattice["geometrically_indecomposable"] == true {
            "conditional_existence"
        } else {
            "unresolved"
        };
        let mut record = json!({"schema_version":"jacobian-certificate/v1",
            "provenance":{"input_blake3":digest(input),"registry_blake3":registry.map(digest),
                "checker_source_blake3":digest(include_bytes!("jacobian_certificate.rs")),
                "algorithm":"maximal-order-minimum-gap/v1","work_budget":budget},
            "candidate":c,"binding":binding(&c,reg.as_ref())?,"lattice":lattice,
            "status":status,"exists_over_base_field":null,
            "required_hypotheses":["E is an ordinary elliptic curve over the declared finite field",
                "the actual geometric endomorphism ring is the specified maximal order",
                "the specified order acts on E over the declared base field; Rosati is conjugation"],
            "hypotheses_verification":"not_performed",
            "curve_equation":null,"cover_maps":null,"dlp_advantage":null});
        let id = format!(
            "urn:jacobian-certificate:1:blake3:{}",
            digest(record.to_string().as_bytes())
        );
        record["certificate_uid"] = json!(id);
        records.push(record);
    }
    Ok(json!({"schema_version":"jacobian-certificates/v1",
        "input_blake3":digest(input),"registry_blake3":registry.map(digest),
        "checker_source_blake3":digest(include_bytes!("jacobian_certificate.rs")),
        "algorithm":"maximal-order-minimum-gap/v1","work_budget":budget,"records":records}))
}

fn sql_literal(value: &str) -> String {
    format!("'{}'", value.replace('\'', "''"))
}
fn optional_sql(v: &Value) -> String {
    v.as_str().map(sql_literal).unwrap_or_else(|| "NULL".into())
}
fn sql(report: &Value) -> String {
    let mut s = String::from("BEGIN;\nCREATE TABLE IF NOT EXISTS jacobian_certificates (certificate_uid TEXT PRIMARY KEY, curve_uid TEXT, curve_slug TEXT, discriminant INTEGER NOT NULL, status TEXT NOT NULL, record_json TEXT NOT NULL CHECK(json_valid(record_json)));\nCREATE INDEX IF NOT EXISTS jacobian_certificates_curve ON jacobian_certificates(curve_uid,status);\nCREATE INDEX IF NOT EXISTS jacobian_certificates_order ON jacobian_certificates(discriminant,status);\n");
    for record in report["records"].as_array().expect("generated records") {
        s.push_str(&format!("INSERT INTO jacobian_certificates VALUES ({},{},{},{},{},{}) ON CONFLICT(certificate_uid) DO NOTHING;\n",
            optional_sql(&record["certificate_uid"]),optional_sql(&record["binding"]["curve_uid"]),
            optional_sql(&record["binding"]["slug"]),record["candidate"]["discriminant"],
            optional_sql(&record["status"]),sql_literal(&record.to_string())));
    }
    s.push_str("COMMIT;\n");
    s
}

#[derive(Parser)]
#[command(about = "Certify bounded CM Hermitian lattices and export conditional Jacobian findings")]
struct Args {
    #[arg(long, default_value = "docs/curves/jacobian-candidates.json")]
    input: PathBuf,
    #[arg(long)]
    registry: Option<PathBuf>,
    #[arg(long, default_value = "docs/curves/jacobian-certificates.json")]
    output: PathBuf,
    /// Optional SQLite migration/import script; no database connection is opened.
    #[arg(long)]
    sql: Option<PathBuf>,
    /// Replay mathematics and require byte-identical JSON (and SQL when supplied).
    #[arg(long)]
    check: bool,
    #[arg(long, default_value_t = 2_000_000)]
    budget: u64,
}

fn same_path(a: &PathBuf, b: &PathBuf) -> bool {
    a == b || matches!((a.canonicalize(),b.canonicalize()),(Ok(x),Ok(y)) if x==y)
}
fn run() -> Check<()> {
    let a = Args::parse();
    let mut paths = vec![&a.input, &a.output];
    if let Some(p) = &a.registry {
        paths.push(p);
    }
    if let Some(p) = &a.sql {
        paths.push(p);
    }
    for i in 0..paths.len() {
        for j in 0..i {
            if same_path(paths[i], paths[j]) {
                return Err("input, registry and outputs must be distinct paths".into());
            }
        }
    }
    let input = fs::read(&a.input).map_err(|e| e.to_string())?;
    let registry = a
        .registry
        .as_ref()
        .map(fs::read)
        .transpose()
        .map_err(|e| e.to_string())?;
    let result = report(&input, registry.as_deref(), a.budget)?;
    let mut outputs = vec![(
        a.output,
        serde_json::to_string_pretty(&result).map_err(|e| e.to_string())? + "\n",
    )];
    if let Some(path) = a.sql {
        outputs.push((path, sql(&result)));
    }
    for (path, text) in outputs {
        if a.check {
            if fs::read_to_string(&path).map_err(|e| e.to_string())? != text {
                return Err(format!("stale or altered certificate: {}", path.display()));
            }
        } else {
            fs::write(path, text).map_err(|e| e.to_string())?;
        }
    }
    println!(
        "{} records; geometric hypotheses remain unverified",
        result["records"].as_array().unwrap().len()
    );
    Ok(())
}
fn main() {
    if let Err(e) = run() {
        eprintln!("jacobian_certificate: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
#[path = "jacobian_certificate/tests.rs"]
mod tests;

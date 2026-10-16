//! Exhaust a finite degree graph without confusing a budget stop with completion.
//! Ordinary F_p classes, 3 < p < 2^61; D_K < -4, |D_K| <= 10^7.
use isogeny_algos::{
    api::{self, PrimeCurve, Status},
    curve::{self, Curve, Pt, RatIsogeny},
    field::{Field, Rng, Zp},
    find::divpoly::kernel_polys,
    icv1,
    int::Int,
    kernel::kohel::kohel,
    poly,
};
use rusqlite::{params, Connection, OptionalExtension};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, HashMap},
    fs,
    path::Path,
};

type Result<T> = std::result::Result<T, String>;
const SCHEMA: &str = "ordinary-isogeny-class/v1";
const HELP: &str = "isogeny-class: finite-degree, complete-vertex isogeny classes
  plan   --curve curve.json [--degrees auto|2,3,5]
  build  --curve curve.json --out fresh.sqlite [--degrees auto|2,3,5] [--steps all|10000]
  resume --input graph.sqlite [--steps all|10000]
  status --input graph.sqlite
  vertex --input graph.sqlite --id ID       exact models, kernels and rational maps
  verify --input graph.sqlite
  solve  --input graph.sqlite --vertex ID --px X --py Y --qx X --qy Y --order R [--max-baby 1000000]
curve.json: {\"p\":113,\"a\":A,\"b\":B,\"order\":112}; order is optional (SEA counts it).
All degrees must be primes <=31 different from p. Auto includes conductor primes
and a Minkowski generating set on the surface. Prime fields below 2^61 only;
ordinary D_K<-4, |D_K|<=10^7. Other inputs fail explicitly.
Every expansion is one SQLite transaction. Resume retains all prior expansions.
Default 10000 expansions is a PARTIAL checkpoint, never a completeness claim.
Use --steps all to exhaust the frontier. Exit 3 means partial; 0 means complete
for build/resume/verify, and a verified logarithm for solve. Plan/status exit 0
reports inspection success, not graph completion. No timings are benchmark claims.
";

fn err<E: std::fmt::Display>(e: E) -> String {
    e.to_string()
}
fn hash(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
fn mul(a: u64, b: u64, n: u64) -> u64 {
    (a as u128 * b as u128 % n as u128) as u64
}
fn pow(mut a: u64, mut e: u64, n: u64) -> u64 {
    let mut r = 1;
    while e > 0 {
        if e & 1 != 0 {
            r = mul(r, a, n);
        }
        a = mul(a, a, n);
        e >>= 1;
    }
    r
}
fn prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'base: for a in [2, 325, 9375, 28178, 450775, 9780504, 1795265022] {
        if a % n == 0 {
            continue;
        }
        let mut x = pow(a % n, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul(x, x, n);
            if x == n - 1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}
fn factor(n: u64) -> Result<Vec<(u64, u32)>> {
    fn split(n: u64, out: &mut Vec<u64>) -> Result<()> {
        if n == 1 {
            return Ok(());
        }
        if prime(n) {
            out.push(n);
            return Ok(());
        }
        if n % 2 == 0 {
            out.push(2);
            return split(n / 2, out);
        }
        for c in 1..=64 {
            let mut x = 2;
            let mut y = 2;
            for _ in 0..500_000 {
                x = ((mul(x, x, n) as u128 + c as u128) % n as u128) as u64;
                for _ in 0..2 {
                    y = ((mul(y, y, n) as u128 + c as u128) % n as u128) as u64;
                }
                let d = gcd(x.abs_diff(y), n);
                if d == n {
                    break;
                }
                if d > 1 {
                    split(d, out)?;
                    return split(n / d, out);
                }
            }
        }
        Err(format!(
            "factorization budget exhausted on {n}; no completeness claim"
        ))
    }
    let mut raw = vec![];
    split(n, &mut raw)?;
    raw.sort_unstable();
    let mut out: Vec<(u64, u32)> = vec![];
    for p in raw {
        if let Some((q, e)) = out.last_mut() {
            if *q == p {
                *e += 1;
                continue;
            }
        }
        out.push((p, 1));
    }
    Ok(out)
}
fn isqrt(n: u64) -> u64 {
    let mut lo = 0;
    let mut hi = 1u64 << 32;
    while lo + 1 < hi {
        let m = (lo + hi) / 2;
        if m as u128 * m as u128 <= n as u128 {
            lo = m;
        } else {
            hi = m;
        }
    }
    lo
}
fn character(d: i64, l: u64) -> i64 {
    if l == 2 {
        return match d.rem_euclid(8) {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => 0,
        };
    }
    let r = (d as i128).rem_euclid(l as i128) as u64;
    if r == 0 {
        0
    } else if pow(r, (l - 1) / 2, l) == 1 {
        1
    } else {
        -1
    }
}
fn fundamental(delta: i64) -> Result<(i64, u64)> {
    let mut square = 1u64;
    let mut sf = 1u64;
    for (p, e) in factor(delta.unsigned_abs())? {
        square *= p.pow(e / 2);
        if e % 2 == 1 {
            sf *= p;
        }
    }
    let d = -(sf as i64);
    if d.rem_euclid(4) == 1 {
        Ok((d, square))
    } else if square % 2 == 0 {
        Ok((4 * d, square / 2))
    } else {
        Err("invalid quadratic discriminant".into())
    }
}
fn class_number(d: i64) -> Result<u64> {
    if d >= -4 || d.unsigned_abs() > 10_000_000 {
        return Err("class-number backend requires -10000000 <= D_K < -4 (special j and larger fields need another backend)".into());
    }
    let limit = isqrt(d.unsigned_abs() / 3) as i64;
    let mut h = 0;
    for a in 1..=limit {
        for b in -a..=a {
            let num = b * b - d;
            if num % (4 * a) != 0 {
                continue;
            }
            let c = num / (4 * a);
            if c < a || ((b.abs() == a || a == c) && b < 0) {
                continue;
            }
            if gcd(gcd(a as u64, b.unsigned_abs()), c as u64) == 1 {
                h += 1;
            }
        }
    }
    Ok(h)
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub p: u64,
    pub a: u64,
    pub b: u64,
    #[serde(default)]
    pub order: Option<u64>,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Stratum {
    pub conductor: u64,
    pub expected_vertices: u64,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Plan {
    pub schema: String,
    pub input: Input,
    pub order: u64,
    pub trace: i64,
    pub discriminant: i64,
    pub fundamental_discriminant: i64,
    pub frobenius_conductor: u64,
    pub conductor_factors: Vec<(u64, u32)>,
    pub maximal_class_number: u64,
    pub strata: Vec<Stratum>,
    pub expected_vertices: u64,
    pub degrees: Vec<u64>,
    pub auto_degrees: Vec<u64>,
    pub order_evidence: Value,
    pub source_identity: Value,
    pub source_code_sha256: String,
    pub storage: Value,
}

fn identities(input: &Input, n: u64) -> Result<Value> {
    let id = icv1::prime(
        &Int::from(input.p),
        &Int::from(input.a),
        &Int::from(input.b),
        &Int::from(n),
    )
    .ok_or("identity failed")?;
    Ok(
        json!({"icv1":id.icv1,"slug":id.slug,"model_json":id.model_json,
        "ec1":null,"curve_uid":null,"ec1_status":"no subgroup/generator selected for this graph vertex; not a cross-method comparison"}),
    )
}
pub fn plan(mut input: Input, degree_text: &str) -> Result<Plan> {
    if input.p <= 3 || input.p >= 1 << 61 || !prime(input.p) {
        return Err("requires a prime 3 < p < 2^61".into());
    }
    input.a %= input.p;
    input.b %= input.p;
    if let Some(n) = input.order {
        let t = (input.p as i128 + 1 - n as i128).unsigned_abs();
        if n == 0 || t * t > 4 * input.p as u128 {
            return Err("supplied group order violates the Hasse interval".into());
        }
    }
    let pc = PrimeCurve::new(Int::from(input.p), Int::from(input.a), Int::from(input.b))?;
    let count = if let Some(n) = input.order {
        api::check_order(&pc, &Int::from(n), 1)
    } else {
        api::count_points(&pc, 1)?
    };
    let n: u64 = count.order.to_string().parse().map_err(err)?;
    let exhaustive_order = if input.p <= 1 << 20 {
        Some(curve::order(
            &Zp::new(input.p),
            &Curve::new(input.a, input.b),
            &mut Rng::new(1),
        ))
    } else {
        None
    };
    if exhaustive_order.is_some_and(|exact| exact != n) {
        return Err("supplied order disagrees with exact point enumeration".into());
    }
    if (exhaustive_order.is_none() && !count.certificate.certified)
        || count
            .certificate
            .checks
            .iter()
            .any(|c| c.status == Status::Fail)
    {
        return Err(
            "group order not certified; point annihilation alone cannot certify the isogeny class"
                .into(),
        );
    }
    let trace = (input.p + 1) as i64 - n as i64;
    if trace == 0 {
        return Err("supersingular input: this backend is ordinary only".into());
    }
    let disc = (trace as i128 * trace as i128 - 4 * input.p as i128) as i64;
    let (dk, f) = fundamental(disc)?;
    let hd = class_number(dk)?;
    let fac = factor(f)?;
    let mut strata = vec![Stratum {
        conductor: 1,
        expected_vertices: hd,
    }];
    for &(l, e) in &fac {
        let previous = strata.clone();
        let mut power = 1;
        for k in 1..=e {
            power *= l;
            let hk = (l as i64 - character(dk, l)) as u64 * l.pow(k - 1);
            for s in &previous {
                strata.push(Stratum {
                    conductor: s.conductor * power,
                    expected_vertices: s.expected_vertices * hk,
                });
            }
        }
    }
    strata.sort_by_key(|s| s.conductor);
    let total = strata.iter().map(|s| s.expected_vertices).sum::<u64>();
    // Every ideal class has an integral ideal of norm <= (2/pi)*sqrt(|D_K|).
    // The larger integer sqrt(|D_K|) is a conservative bound, no floating point.
    let mut auto = fac.iter().map(|x| x.0).collect::<Vec<_>>();
    if hd > 1 {
        for l in 2..=isqrt(dk.unsigned_abs()) {
            if prime(l) && character(dk, l) != -1 {
                auto.push(l);
            }
        }
    }
    auto.sort_unstable();
    auto.dedup();
    let mut degrees = if degree_text == "auto" {
        auto.clone()
    } else if degree_text.is_empty() {
        vec![]
    } else {
        degree_text
            .split(',')
            .map(|x| x.parse::<u64>().map_err(err))
            .collect::<Result<Vec<_>>>()?
    };
    degrees.sort_unstable();
    degrees.dedup();
    if degrees.iter().any(|&l| !prime(l) || l == input.p || l > 31) {
        return Err(format!(
            "degree backend requires primes <=31, different from p; automatic set is {auto:?}"
        ));
    }
    let evidence = json!({"certified":true,"exact_small_field_enumeration":exhaustive_order,"via_twist":count.certificate.via_twist,
        "divisor":count.certificate.divisor.to_string(),"twist_divisor":count.twist_certificate.divisor.to_string(),
        "order_method":format!("{:?}",count.method),"point_samples":count.certificate.points.iter().map(|(x,y)|json!([x.to_string(),y.to_string()])).collect::<Vec<_>>(),
        "factor_primality":"native deterministic u64 Miller-Rabin rechecked for discriminant; order certificate uses isogeny_algos certificate backend"});
    let identity = identities(&input, n)?;
    Ok(Plan {
        schema: SCHEMA.into(),
        input,
        order: n,
        trace,
        discriminant: disc,
        fundamental_discriminant: dk,
        frobenius_conductor: f,
        conductor_factors: fac,
        maximal_class_number: hd,
        strata,
        expected_vertices: total,
        degrees,
        auto_degrees: auto,
        order_evidence: evidence,
        source_identity: identity,
        source_code_sha256: hash(include_bytes!("lib.rs")),
        storage: json!({"fixed_u64_vertex_payload_bytes":total as u128*24,
            "layout_includes":"a,b,j only (3 u64); excludes indexes, maps, journals, free space and DLP evidence; compressed or generative encodings can be smaller",
            "sqlite_bytes":"not predicted; depends on degree set and SQLite layout","bounded_ram":"SQLite pages and one vertex's outgoing kernels; no whole-class HashMap"}),
    })
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
struct Edge {
    target: u64,
    ell: u64,
    kernel: Vec<u64>,
    u: u64,
}
fn kernels(f: &Zp, c: &Curve<u64>, ell: u64) -> Vec<Vec<u64>> {
    let mut rng = Rng::new(1);
    let mut hs = if ell == 2 {
        poly::roots(f, &vec![c.b, c.a, 0, 1], &mut rng)
            .into_iter()
            .map(|r| vec![f.neg(r), 1])
            .collect()
    } else {
        kernel_polys(f, c, ell, &mut rng)
    };
    hs.sort();
    hs
}
fn iso_scale(f: &Zp, raw: &Curve<u64>, target: &Curve<u64>) -> Result<u64> {
    if raw.a == 0 || raw.b == 0 || target.a == 0 || target.b == 0 {
        return Err("special j model encountered".into());
    }
    let u2 = f.div(f.mul(target.b, raw.a), f.mul(raw.b, target.a));
    let u = f
        .sqrt(u2)
        .ok_or("same-j codomain is a different F_p twist")?;
    let u = u.min(f.neg(u));
    if f.mul(f.sq(u2), raw.a) != target.a || f.mul(f.mul(f.sq(u2), u2), raw.b) != target.b {
        return Err("codomain isomorphism failed".into());
    }
    Ok(u)
}
fn map_identity(f: &Zp, m: &RatIsogeny<Zp>) -> bool {
    let d = &m.den;
    let n = &m.num;
    let deriv = poly::sub(
        f,
        &poly::mul(f, &poly::derivative(f, n), d),
        &poly::mul(f, n, &poly::derivative(f, d)),
    );
    let lhs = poly::mul(
        f,
        &vec![m.dom.b, m.dom.a, 0, 1],
        &poly::mul(f, &deriv, &deriv),
    );
    let d2 = poly::mul(f, d, d);
    let d3 = poly::mul(f, &d2, d);
    let d4 = poly::mul(f, &d2, &d2);
    let n3 = poly::mul(f, &poly::mul(f, n, n), n);
    let rhs = poly::add(
        f,
        &poly::add(
            f,
            &poly::mul(f, &n3, d),
            &poly::scale(f, &poly::mul(f, n, &d3), m.cod.a),
        ),
        &poly::scale(f, &d4, m.cod.b),
    );
    poly::sub(f, &lhs, &rhs).is_empty()
}
fn open(path: &Path) -> Result<Connection> {
    if !path.is_file() {
        return Err("input database does not exist".into());
    }
    let db = Connection::open(path).map_err(err)?;
    db.busy_timeout(std::time::Duration::from_secs(30))
        .map_err(err)?;
    db.execute_batch("PRAGMA foreign_keys=ON; PRAGMA synchronous=FULL; PRAGMA cache_size=-32768;")
        .map_err(err)?;
    Ok(db)
}
fn load(db: &Connection) -> Result<Plan> {
    let s: String = db
        .query_row("SELECT value FROM metadata WHERE key='plan'", [], |r| {
            r.get(0)
        })
        .map_err(err)?;
    let plan: Plan = serde_json::from_str(&s).map_err(err)?;
    if plan.schema != SCHEMA {
        return Err("unsupported database schema".into());
    }
    Ok(plan)
}
fn counts(db: &Connection) -> Result<(u64, u64, u64)> {
    Ok((
        db.query_row("SELECT count(*) FROM nodes", [], |r| r.get(0))
            .map_err(err)?,
        db.query_row("SELECT count(*) FROM nodes WHERE expanded=1", [], |r| {
            r.get(0)
        })
        .map_err(err)?,
        db.query_row("SELECT count(*) FROM edges", [], |r| r.get(0))
            .map_err(err)?,
    ))
}
fn create(path: &Path, p: &Plan) -> Result<Connection> {
    // create_new prevents a race from overwriting an existing database.
    fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(err)?;
    let db = open(path)?;
    db.execute_batch("CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE nodes(id INTEGER PRIMARY KEY,a INTEGER NOT NULL,b INTEGER NOT NULL,j INTEGER NOT NULL UNIQUE,expanded INTEGER NOT NULL DEFAULT 0,conductor INTEGER);
        CREATE INDEX frontier ON nodes(expanded,id);
        CREATE TABLE edges(source INTEGER NOT NULL REFERENCES nodes(id),ell INTEGER NOT NULL,kernel TEXT NOT NULL,target INTEGER NOT NULL REFERENCES nodes(id),u INTEGER NOT NULL,PRIMARY KEY(source,ell,kernel)) WITHOUT ROWID;
        CREATE INDEX reverse_edges ON edges(ell,target,source);").map_err(err)?;
    db.execute(
        "INSERT INTO metadata VALUES('plan',?)",
        [serde_json::to_string(p).map_err(err)?],
    )
    .map_err(err)?;
    let f = Zp::new(p.input.p);
    let c = Curve::new(p.input.a, p.input.b);
    db.execute(
        "INSERT INTO nodes(id,a,b,j) VALUES(0,?,?,?)",
        params![c.a, c.b, curve::jinv(&f, &c)],
    )
    .map_err(err)?;
    Ok(db)
}

fn expand(db: &mut Connection, p: &Plan, steps: u64) -> Result<()> {
    if p.source_code_sha256 != hash(include_bytes!("lib.rs")) {
        return Err(
            "source code differs from frozen plan; replay/migrate explicitly before resuming"
                .into(),
        );
    }
    let f = Zp::new(p.input.p);
    for _ in 0..steps {
        let next: Option<(u64, u64, u64)> = db
            .query_row(
                "SELECT id,a,b FROM nodes WHERE expanded=0 ORDER BY id LIMIT 1",
                [],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
            )
            .optional()
            .map_err(err)?;
        let Some((id, a, b)) = next else {
            break;
        };
        let c = Curve::new(a, b);
        let tx = db.transaction().map_err(err)?;
        for &ell in &p.degrees {
            for h in kernels(&f, &c, ell) {
                let m = kohel(&f, &c, &h, ell);
                if !map_identity(&f, &m) {
                    return Err(format!("map identity failed at vertex {id}, degree {ell}"));
                }
                let j = curve::jinv(&f, &m.cod);
                tx.execute(
                    "INSERT OR IGNORE INTO nodes(a,b,j) VALUES(?,?,?)",
                    params![m.cod.a, m.cod.b, j],
                )
                .map_err(err)?;
                let (target, ta, tb): (u64, u64, u64) = tx
                    .query_row("SELECT id,a,b FROM nodes WHERE j=?", [j], |r| {
                        Ok((r.get(0)?, r.get(1)?, r.get(2)?))
                    })
                    .map_err(err)?;
                let u = iso_scale(&f, &m.cod, &Curve::new(ta, tb))?;
                tx.execute(
                    "INSERT INTO edges VALUES(?,?,?,?,?)",
                    params![id, ell, serde_json::to_string(&h).map_err(err)?, target, u],
                )
                .map_err(err)?;
            }
        }
        tx.execute("UPDATE nodes SET expanded=1 WHERE id=?", [id])
            .map_err(err)?;
        tx.commit().map_err(err)?;
    }
    Ok(())
}

fn stratify(db: &mut Connection, p: &Plan) -> Result<Value> {
    let (n, expanded, _) = counts(db)?;
    if expanded != n || n != p.expected_vertices {
        return Ok(
            json!({"status":"NOT_ESTABLISHED","reason":"requires exhausted frontier and exact class population"}),
        );
    }
    if p.conductor_factors
        .iter()
        .any(|(ell, _)| !p.degrees.contains(ell))
    {
        return Ok(json!({"status":"NOT_ESTABLISHED","reason":"missing conductor prime degree"}));
    }
    let tx = db.transaction().map_err(err)?;
    tx.execute("UPDATE nodes SET conductor=1", [])
        .map_err(err)?;
    for &(ell, depth) in &p.conductor_factors {
        tx.execute_batch("DROP TABLE IF EXISTS temp.heights; CREATE TEMP TABLE heights(id INTEGER PRIMARY KEY,distance INTEGER NOT NULL);").map_err(err)?;
        tx.execute("INSERT INTO heights SELECT source,0 FROM edges WHERE ell=? GROUP BY source HAVING count(*)=1",[ell]).map_err(err)?;
        for k in 0..depth {
            tx.execute("INSERT OR IGNORE INTO heights SELECT e.target,? FROM heights h JOIN edges e ON e.source=h.id WHERE h.distance=? AND e.ell=?",params![k+1,k,ell]).map_err(err)?;
        }
        let reached: u64 = tx
            .query_row("SELECT count(*) FROM heights", [], |r| r.get(0))
            .map_err(err)?;
        if reached != n {
            return Err(format!(
                "degree {ell} floor-distance labeling reached {reached}/{n}"
            ));
        }
        for distance in 0..=depth {
            tx.execute("UPDATE nodes SET conductor=conductor*? WHERE id IN (SELECT id FROM heights WHERE distance=?)",params![ell.pow(depth-distance),distance]).map_err(err)?;
        }
    }
    let observed = {
        let mut stmt = tx
            .prepare("SELECT conductor,count(*) FROM nodes GROUP BY conductor ORDER BY conductor")
            .map_err(err)?;
        let rows = stmt
            .query_map([], |r| {
                Ok(Stratum {
                    conductor: r.get(0)?,
                    expected_vertices: r.get(1)?,
                })
            })
            .map_err(err)?;
        rows.collect::<std::result::Result<Vec<_>, _>>()
            .map_err(err)?
    };
    if observed != p.strata {
        return Err("stratum populations disagree with exact ring-class counts".into());
    }
    tx.commit().map_err(err)?;
    Ok(
        json!({"status":"VERIFIED_POPULATIONS","strata":observed,"method":"prime-degree floor distance; ordinary volcano theorem, D_K<-4"}),
    )
}
fn summary(db: &Connection, p: &Plan) -> Result<Value> {
    let (n, expanded, e) = counts(db)?;
    let graph_complete = n == expanded;
    let class_complete = graph_complete && n == p.expected_vertices;
    let mut statement = db
        .prepare("SELECT ell,count(*) FROM edges GROUP BY ell ORDER BY ell")
        .map_err(err)?;
    let degrees = statement
        .query_map([], |r| Ok((r.get::<_, u64>(0)?, r.get::<_, u64>(1)?)))
        .map_err(err)?
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(err)?;
    let validated: Option<String> = db
        .query_row(
            "SELECT value FROM metadata WHERE key='verification'",
            [],
            |r| r.get(0),
        )
        .optional()
        .map_err(err)?;
    Ok(
        json!({"schema":SCHEMA,"status":if class_complete {"COMPLETE_SELECTED_DEGREE_CLASS"} else if graph_complete {"COMPLETE_COMPONENT_ONLY"} else {"PARTIAL"},
        "vertices":n,"expanded_vertices":expanded,"frontier_vertices":n-expanded,"expected_class_vertices":p.expected_vertices,
        "directed_kernel_edges":e,"edges_by_degree":degrees,"selected_degrees":p.degrees,
        "vertex_class_complete":class_complete,"selected_degree_edges_complete":graph_complete,
        "all_arbitrary_degree_maps":"INFINITE_FAMILY_NOT_AN_ENUMERATION_TARGET",
        "independent_implementation_replay":"NOT_RUN","replay_receipt":validated.and_then(|s|serde_json::from_str::<Value>(&s).ok()),
        "dlp_workload":"NOT_RUN_BY_GRAPH_BUILD","expected_strata":p.strata,"source_identity":p.source_identity}),
    )
}

fn read_edges(db: &Connection, id: u64) -> Result<Vec<Edge>> {
    let mut s = db
        .prepare("SELECT target,ell,kernel,u FROM edges WHERE source=? ORDER BY ell,kernel")
        .map_err(err)?;
    let raw = s
        .query_map([id], |r| {
            Ok((
                r.get::<_, u64>(0)?,
                r.get::<_, u64>(1)?,
                r.get::<_, String>(2)?,
                r.get::<_, u64>(3)?,
            ))
        })
        .map_err(err)?;
    raw.map(|row| {
        let (t, l, h, u) = row.map_err(err)?;
        Ok(Edge {
            target: t,
            ell: l,
            kernel: serde_json::from_str(&h).map_err(err)?,
            u,
        })
    })
    .collect()
}
fn node(db: &Connection, id: u64) -> Result<Curve<u64>> {
    db.query_row("SELECT a,b FROM nodes WHERE id=?", [id], |r| {
        Ok(Curve::new(r.get(0)?, r.get(1)?))
    })
    .map_err(err)
}

fn describe_vertex(db: &Connection, p: &Plan, id: u64) -> Result<Value> {
    let f = Zp::new(p.input.p);
    let c = node(db, id)?;
    let conductor: Option<u64> = db
        .query_row("SELECT conductor FROM nodes WHERE id=?", [id], |r| r.get(0))
        .map_err(err)?;
    let input = Input {
        p: p.input.p,
        a: c.a,
        b: c.b,
        order: Some(p.order),
    };
    let mut maps = vec![];
    for edge in read_edges(db, id)? {
        let m = kohel(&f, &c, &edge.kernel, edge.ell);
        let target = node(db, edge.target)?;
        let target_input = Input {
            p: p.input.p,
            a: target.a,
            b: target.b,
            order: Some(p.order),
        };
        maps.push(json!({"target_vertex":edge.target,"degree":edge.ell,"kernel_polynomial_ascending":edge.kernel,
            "normalized_quotient":{"a":m.cod.a,"b":m.cod.b},"x_numerator_ascending":m.num,"x_denominator_ascending":m.den,
            "target_model":{"a":target.a,"b":target.b},"target_identity":identities(&target_input,p.order)?,
            "target_isomorphism_u":edge.u,"formula":"X=u^2*N(x)/D(x); Y=u^3*y*(N'(x)D(x)-N(x)D'(x))/D(x)^2",
            "exceptional_inputs":"infinity and kernel map to infinity","field_p":p.input.p}));
    }
    Ok(
        json!({"vertex":id,"model":input,"j":curve::jinv(&f,&c),"endomorphism_conductor":conductor,
        "identity":identities(&input,p.order)?,"outgoing_isogenies":maps,
        "status":"STORED_RECORD; use verify for whole-database replay"}),
    )
}
fn replay(db: &mut Connection, p: &Plan) -> Result<Value> {
    let integrity: String = db
        .query_row("PRAGMA integrity_check", [], |r| r.get(0))
        .map_err(err)?;
    if integrity != "ok" {
        return Err(format!("SQLite integrity_check: {integrity}"));
    }
    let mut fk = db.prepare("PRAGMA foreign_key_check").map_err(err)?;
    if fk.query([]).map_err(err)?.next().map_err(err)?.is_some() {
        return Err("foreign key integrity failed".into());
    }
    drop(fk);
    let recomputed = plan(
        p.input.clone(),
        &p.degrees
            .iter()
            .map(u64::to_string)
            .collect::<Vec<_>>()
            .join(","),
    )?;
    // Source hash may differ for a replay, but all mathematical/input fields must agree.
    let mut saved = p.clone();
    saved.source_code_sha256 = recomputed.source_code_sha256.clone();
    if saved != recomputed {
        return Err("stored plan does not match recomputed class invariants".into());
    }
    let f = Zp::new(p.input.p);
    let mut verified = 0u64;
    let mut id = -1i64;
    loop {
        let row: Option<(u64, u64, u64, u64, u64)> = db
            .query_row(
                "SELECT id,a,b,j,expanded FROM nodes WHERE id>? ORDER BY id LIMIT 1",
                [id],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?, r.get(4)?)),
            )
            .optional()
            .map_err(err)?;
        let Some((index, a, b, j, expanded)) = row else {
            break;
        };
        id = index as i64;
        if a >= f.p || b >= f.p {
            return Err("noncanonical coefficient".into());
        }
        let c = Curve::new(a, b);
        if !curve::is_smooth(&f, &c) || curve::jinv(&f, &c) != j {
            return Err(format!("invalid vertex {index}"));
        }
        if index == 0 && (a != p.input.a || b != p.input.b) {
            return Err("root model mismatch".into());
        }
        if index != 0 {
            let incoming: u64 = db
                .query_row(
                    "SELECT count(*) FROM edges WHERE target=? AND source<?",
                    params![index, index],
                    |r| r.get(0),
                )
                .map_err(err)?;
            if incoming == 0 {
                return Err(format!("vertex {index} has no earlier certified parent"));
            }
        }
        let edges = read_edges(db, index)?;
        if expanded == 0 {
            if !edges.is_empty() {
                return Err("unexpanded vertex has edges".into());
            }
            continue;
        }
        if expanded != 1 {
            return Err("invalid expanded flag".into());
        }
        let mut expected = BTreeMap::new();
        for &ell in &p.degrees {
            for h in kernels(&f, &c, ell) {
                expected.insert((ell, h), false);
            }
        }
        for edge in edges {
            let seen = expected
                .get_mut(&(edge.ell, edge.kernel.clone()))
                .ok_or("unexpected outgoing kernel")?;
            if *seen {
                return Err("duplicate kernel".into());
            }
            *seen = true;
            let m = kohel(&f, &c, &edge.kernel, edge.ell);
            if !map_identity(&f, &m) {
                return Err("exact rational-map polynomial identity failed".into());
            }
            let target = node(db, edge.target)?;
            if edge.u >= f.p
                || f.mul(f.pow(edge.u, 4), m.cod.a) != target.a
                || f.mul(f.pow(edge.u, 6), m.cod.b) != target.b
                || edge.u == 0
            {
                return Err("codomain isomorphism replay failed".into());
            }
            verified += 1;
        }
        if expected.values().any(|seen| !*seen) {
            return Err(format!("missing outgoing kernel at vertex {index}"));
        }
    }
    let strata = stratify(db, p)?;
    let (n, expanded, e) = counts(db)?;
    if verified != e {
        return Err("edge count mismatch".into());
    }
    let receipt = json!({"status":"REPLAYED","vertices":n,"expanded_vertices":expanded,"exact_map_identities":verified,
        "kernel_sets":"regenerated for every expanded vertex","strata":strata,
        "implementation_independence":"shared isogeny_algos arithmetic; independent implementation audit is a separate obligation",
        "graph_complete":n==expanded,"vertex_class_complete":n==expanded&&n==p.expected_vertices,
        "replay_source_sha256":hash(include_bytes!("lib.rs"))});
    db.execute(
        "INSERT OR REPLACE INTO metadata VALUES('verification',?)",
        [serde_json::to_string(&receipt).map_err(err)?],
    )
    .map_err(err)?;
    Ok(receipt)
}

struct Options(BTreeMap<String, String>);
impl Options {
    fn parse(args: &[String], allowed: &[&str]) -> Result<Self> {
        if args.len() % 2 != 0 {
            return Err("every option requires a value".into());
        }
        let mut map = BTreeMap::new();
        for pair in args.chunks(2) {
            let name = pair[0]
                .strip_prefix("--")
                .ok_or("expected --option value")?;
            if !allowed.contains(&name) {
                return Err(format!("unknown option --{name}"));
            }
            if map.insert(name.into(), pair[1].clone()).is_some() {
                return Err(format!("duplicate option --{name}"));
            }
        }
        Ok(Self(map))
    }
    fn get(&self, k: &str) -> Option<&str> {
        self.0.get(k).map(String::as_str)
    }
    fn need(&self, k: &str) -> Result<&str> {
        self.get(k).ok_or_else(|| format!("--{k} is required"))
    }
    fn number(&self, k: &str) -> Result<u64> {
        self.need(k)?.parse().map_err(err)
    }
    fn steps(&self) -> Result<u64> {
        match self.get("steps").unwrap_or("10000") {
            "all" => Ok(u64::MAX),
            s => s.parse().map_err(err),
        }
    }
}

fn solve(db: &Connection, p: &Plan, o: &Options) -> Result<Value> {
    let f = Zp::new(p.input.p);
    let id = o.number("vertex")?;
    let c = node(db, id)?;
    let coordinates = [
        o.number("px")?,
        o.number("py")?,
        o.number("qx")?,
        o.number("qy")?,
    ];
    if coordinates.iter().any(|&x| x >= f.p) {
        return Err("noncanonical point coordinate".into());
    }
    let g = Pt::Aff(coordinates[0], coordinates[1]);
    let q = Pt::Aff(coordinates[2], coordinates[3]);
    let r = o.number("order")?;
    if r < 2
        || r > p.order
        || p.order % r != 0
        || !curve::on_curve(&f, &c, &g)
        || !curve::on_curve(&f, &c, &q)
    {
        return Err("invalid subgroup order or point".into());
    }
    if curve::pmul(&f, &c, &g, r as u128) != Pt::Inf {
        return Err("generator not killed by order".into());
    }
    for (l, _) in factor(r)? {
        if curve::pmul(&f, &c, &g, (r / l) as u128) == Pt::Inf {
            return Err("supplied generator has smaller order".into());
        }
    }
    let mut m = isqrt(r);
    if m * m < r {
        m += 1;
    }
    let cap = o
        .get("max-baby")
        .unwrap_or("1000000")
        .parse::<u64>()
        .map_err(err)?;
    if m > cap {
        return Err(format!(
            "BSGS needs {m} baby entries, --max-baby allows {cap}"
        ));
    }
    let mut table = HashMap::new();
    let mut cur = Pt::Inf;
    for j in 0..m {
        table.entry(cur).or_insert(j);
        cur = curve::padd(&f, &c, &cur, &g);
    }
    let jump = curve::neg(&f, &curve::pmul(&f, &c, &g, m as u128));
    let mut cur = q;
    let mut lookups = 0;
    for i in 0..=m {
        lookups += 1;
        if let Some(&j) = table.get(&cur) {
            let k = (i as u128 * m as u128 + j as u128) % r as u128;
            if curve::pmul(&f, &c, &g, k) == q {
                return Ok(
                    json!({"status":"VERIFIED_LOGARITHM","method":"BSGS","vertex":id,"order":r,"scalar":k,
                "verified_equation":"[scalar]P=Q","baby_entries":m,"giant_lookups":lookups,
                "dlp_scope":"one explicit cyclic-subgroup target; not an all-vertex workload or benchmark",
                "timing":null,"total_operations":null,"operation_accounting":"table sizes/lookups only; scalar checks, field work and allocation not priced"}),
                );
            }
        }
        cur = curve::padd(&f, &c, &cur, &jump);
    }
    Err("target is not in the specified cyclic subgroup".into())
}

/// Shared by `crypto isogeny class`, `isogeny class`, and `isogeny-class`.
pub fn run(args: &[String]) -> Result<(Value, i32)> {
    let command = args.first().map(String::as_str).unwrap_or("help");
    if matches!(command, "help" | "--help" | "-h") {
        return Ok((json!({"help":HELP}), 0));
    }
    match command {
        "plan" | "build" => {
            let o = Options::parse(
                &args[1..],
                if command == "build" {
                    &["curve", "degrees", "out", "steps"]
                } else {
                    &["curve", "degrees"]
                },
            )?;
            let input = serde_json::from_slice::<Input>(&fs::read(o.need("curve")?).map_err(err)?)
                .map_err(err)?;
            let p = plan(input, o.get("degrees").unwrap_or("auto"))?;
            if command == "plan" {
                return Ok((serde_json::to_value(p).map_err(err)?, 0));
            }
            let steps = o.steps()?;
            let path = Path::new(o.need("out")?);
            let mut db = create(path, &p)?;
            expand(&mut db, &p, steps)?;
            let (n, expanded, _) = counts(&db)?;
            if n == expanded && n == p.expected_vertices {
                replay(&mut db, &p)?;
            }
            let s = summary(&db, &p)?;
            let code = if s["vertex_class_complete"] == true {
                0
            } else {
                3
            };
            Ok((s, code))
        }
        "resume" | "status" | "verify" | "solve" | "vertex" => {
            let allowed: &[&str] = match command {
                "resume" => &["input", "steps"],
                "vertex" => &["input", "id"],
                "solve" => &[
                    "input", "vertex", "px", "py", "qx", "qy", "order", "max-baby",
                ],
                _ => &["input"],
            };
            let o = Options::parse(&args[1..], allowed)?;
            let mut db = open(Path::new(o.need("input")?))?;
            let p = load(&db)?;
            if command == "vertex" {
                return Ok((describe_vertex(&db, &p, o.number("id")?)?, 0));
            }
            if command == "solve" {
                return Ok((solve(&db, &p, &o)?, 0));
            }
            if command == "resume" {
                let steps = o.steps()?;
                db.execute("DELETE FROM metadata WHERE key='verification'", [])
                    .map_err(err)?;
                expand(&mut db, &p, steps)?;
                let (n, expanded, _) = counts(&db)?;
                if n == expanded && n == p.expected_vertices {
                    replay(&mut db, &p)?;
                }
            }
            if command == "verify" {
                replay(&mut db, &p)?;
            }
            let s = summary(&db, &p)?;
            let code = if command == "status" || s["vertex_class_complete"] == true {
                0
            } else {
                3
            };
            Ok((s, code))
        }
        _ => Err(format!("unknown command {command}; use --help")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    static NEXT_TEST: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    fn control(p: u64, j: i64) -> Input {
        let f = Zp::new(p);
        let c = curve::from_j(&f, j.rem_euclid(p as i64) as u64);
        let order = curve::order(&f, &c, &mut Rng::new(1));
        Input {
            p,
            a: c.a,
            b: c.b,
            order: Some(order),
        }
    }
    fn database(input: Input) -> (Connection, Plan) {
        let p = plan(input, "auto").unwrap();
        let path = std::env::temp_dir().join(format!(
            "volcano-control-{}-{}-{}.sqlite",
            std::process::id(),
            p.input.p,
            NEXT_TEST.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        let _ = fs::remove_file(&path);
        let db = create(&path, &p).unwrap();
        (db, p)
    }
    #[test]
    fn complete_controls_match_independent_model_census_and_strata() {
        for (p, j, total) in [(113, -3375, 8), (223, -32768, 9), (67, -3375, 10)] {
            let input = control(p, j);
            let (mut db, plan) = database(input);
            assert_eq!(plan.expected_vertices, total);
            expand(&mut db, &plan, u64::MAX).unwrap();
            let receipt = replay(&mut db, &plan).unwrap();
            assert_eq!(receipt["vertex_class_complete"], true);
            // Independent exhaustive census of every nonsingular model, with
            // literal y-square enumeration rather than an isogeny traversal.
            let mut js = std::collections::BTreeSet::new();
            let f = Zp::new(p);
            let mut multiplicities = vec![0u64; p as usize];
            for y in 0..p {
                multiplicities[mul(y, y, p) as usize] += 1;
            }
            for a in 0..p {
                for b in 0..p {
                    let c = Curve::new(a, b);
                    if !curve::is_smooth(&f, &c) {
                        continue;
                    }
                    let n = 1
                        + (0..p)
                            .map(|x| multiplicities[curve::rhs(&f, &c, x) as usize])
                            .sum::<u64>();
                    if n == plan.order {
                        js.insert(curve::jinv(&f, &c));
                    }
                }
            }
            let saved = db
                .prepare("SELECT j FROM nodes")
                .unwrap()
                .query_map([], |r| r.get(0))
                .unwrap()
                .collect::<std::result::Result<std::collections::BTreeSet<u64>, _>>()
                .unwrap();
            assert_eq!(saved, js);
            assert_eq!(js.len() as u64, total);
        }
    }
    #[test]
    fn resume_equals_uninterrupted_and_tampering_fails() {
        let input = control(67, -3375);
        let (mut db, p) = database(input.clone());
        expand(&mut db, &p, 1).unwrap();
        assert_eq!(summary(&db, &p).unwrap()["vertex_class_complete"], false);
        replay(&mut db, &p).unwrap();
        expand(&mut db, &p, u64::MAX).unwrap();
        replay(&mut db, &p).unwrap();
        let n = counts(&db).unwrap();
        assert_eq!(n.0, 10);
        assert_eq!(n.1, 10);
        let (mut uninterrupted, reference_plan) = database(input);
        expand(&mut uninterrupted, &reference_plan, u64::MAX).unwrap();
        assert_eq!(counts(&uninterrupted).unwrap(), n);
        for id in 0..n.0 {
            assert_eq!(node(&db, id).unwrap(), node(&uninterrupted, id).unwrap());
            assert_eq!(
                read_edges(&db, id).unwrap(),
                read_edges(&uninterrupted, id).unwrap()
            );
        }
        db.execute("UPDATE edges SET u=0 WHERE source=0", [])
            .unwrap();
        assert!(replay(&mut db, &p).is_err());
    }
    #[test]
    fn finite_component_is_not_a_complete_class() {
        let input = control(67, -3375);
        let p = plan(input, "2").unwrap();
        let path =
            std::env::temp_dir().join(format!("volcano-component-{}.sqlite", std::process::id()));
        let _ = fs::remove_file(&path);
        let mut db = create(&path, &p).unwrap();
        expand(&mut db, &p, u64::MAX).unwrap();
        let s = summary(&db, &p).unwrap();
        assert_eq!(s["selected_degree_edges_complete"], true);
        assert_eq!(s["vertex_class_complete"], false);
    }
    #[test]
    fn solve_checks_the_exact_cyclic_subgroup_and_budget() {
        let (db, p) = database(control(113, -3375));
        let f = Zp::new(p.input.p);
        let c = Curve::new(p.input.a, p.input.b);
        let mut rng = Rng::new(9);
        let subgroup = factor(p.order).unwrap().last().unwrap().0;
        let g = loop {
            let g = curve::pmul(
                &f,
                &c,
                &c.random_point(&f, &mut rng),
                (p.order / subgroup) as u128,
            );
            if g != Pt::Inf {
                break g;
            }
        };
        let q = curve::pmul(&f, &c, &g, 5);
        let Pt::Aff(x, y) = g else { unreachable!() };
        let Pt::Aff(qx, qy) = q else { unreachable!() };
        let mut opts = Options(BTreeMap::from([
            ("vertex".into(), "0".into()),
            ("px".into(), x.to_string()),
            ("py".into(), y.to_string()),
            ("qx".into(), qx.to_string()),
            ("qy".into(), qy.to_string()),
            ("order".into(), subgroup.to_string()),
        ]));
        assert_eq!(solve(&db, &p, &opts).unwrap()["scalar"], 5);
        opts.0.insert("max-baby".into(), "1".into());
        assert!(solve(&db, &p, &opts).is_err());
        opts.0.remove("max-baby");
        let r = poly::roots(&f, &vec![c.b, c.a, 0, 1], &mut rng)[0];
        opts.0.insert("qx".into(), r.to_string());
        opts.0.insert("qy".into(), "0".into());
        assert!(solve(&db, &p, &opts).is_err());
    }
    #[test]
    fn automatic_degrees_connect_a_nontrivial_surface_class_group() {
        let f = Zp::new(59);
        let mut input = None;
        'search: for a in 0..59 {
            for b in 0..59 {
                let c = Curve::new(a, b);
                if curve::is_smooth(&f, &c) && curve::order(&f, &c, &mut Rng::new(1)) == 48 {
                    input = Some(Input {
                        p: 59,
                        a,
                        b,
                        order: Some(48),
                    });
                    break 'search;
                }
            }
        }
        let (mut db, p) = database(input.unwrap());
        assert_eq!(p.fundamental_discriminant, -23);
        assert_eq!(p.maximal_class_number, 3);
        assert_eq!(p.degrees, vec![2, 3]);
        assert_eq!(p.expected_vertices, 6);
        expand(&mut db, &p, u64::MAX).unwrap();
        assert_eq!(replay(&mut db, &p).unwrap()["vertex_class_complete"], true);
    }
    #[test]
    fn domain_and_budget_rejections() {
        let mut bad = control(113, -3375);
        bad.order = Some(0);
        assert!(plan(bad, "auto").is_err());
        let mut oversized = control(113, -3375);
        oversized.order = Some(u64::MAX);
        assert!(plan(oversized, "auto").is_err());
        assert!(plan(
            Input {
                p: 9,
                a: 1,
                b: 1,
                order: None
            },
            "auto"
        )
        .is_err());
        assert!(plan(control(113, -3375), "4").is_err());
        assert!(plan(control(113, -3375), "113").is_err());
        for n in [2, 3, 5, 7, 11, 65537, 716481758899867] {
            assert!(prime(n));
        }
        for n in [0, 1, 9, 341, 3215031751] {
            assert!(!prime(n));
        }
        assert_eq!(factor(2u64.pow(28) * 7).unwrap(), vec![(2, 28), (7, 1)]);
    }
}

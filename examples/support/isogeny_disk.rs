//! Transactional BFS: memory contains one source and its neighbours, not the class.
use super::disk_sqlite::{quote, Db};
#[path = "velu_transport.rs"]
mod velu_transport;
#[cfg(test)]
use crypto_lib::cryptanalysis::isogeny_walk::walk::StartCurve;
use crypto_lib::cryptanalysis::{
    isogeny_degree_search::class_number::{factor_u64, is_prime_u64},
    isogeny_walk::{
        curve::{self, Model},
        field::Field,
        kernel,
        modpoly::{self, ModPoly},
        poly,
        record::{ec1_identity, ec1_records},
    },
    prime_fast::{collect_relations, solve_relations, FastCurve, FastFactorBase, FastRelation},
};
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    fs,
    io::{self, Write},
    path::Path,
};
use velu_transport::KernelMap;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Spec {
    p: u64,
    r: u64,
    a: u64,
    b: u64,
    g: [u64; 2],
    q: [u64; 2],
    public_check_scalar: u64,
    primes: Vec<u64>,
    runs: u64,
    seed_start: u64,
    trials: u64,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Node {
    j: u64,
    a: u64,
    b: u64,
    g: [u64; 2],
    q: [u64; 2],
    depth: u64,
    parent_edge: Option<i64>,
    ec1: String,
    uid: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Edge {
    source: i64,
    target: i64,
    degree: u64,
    kernel: Vec<u64>,
    u2: u64,
    u3: u64,
    image_g: [u64; 2],
    image_q: [u64; 2],
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Relation {
    a: u64,
    b: u64,
    entries: Vec<(u32, i8)>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Run {
    seed: u64,
    targets: u64,
    relations: Vec<Relation>,
    rejected: usize,
    scalar: Option<u64>,
    status: String,
}

fn word(f: &Field, a: &crypto_lib::cryptanalysis::isogeny_walk::field::Fe) -> u64 {
    f.to_big(a).to_u64().unwrap()
}
fn model(f: &Field, n: &Node) -> Model {
    Model {
        a: f.from_u64(n.a),
        b: f.from_u64(n.b),
    }
}
fn pair(p: Option<(u64, u64)>) -> [u64; 2] {
    let (x, y) = p.expect("nonzero subgroup image");
    [x, y]
}
fn tuple(p: [u64; 2]) -> Option<(u64, u64)> {
    Some((p[0], p[1]))
}
fn encode<T: Serialize>(v: &T) -> String {
    serde_json::to_string(v).unwrap()
}
fn read_node(db: &Db, id: i64) -> Node {
    serde_json::from_str(&db.text(&format!("SELECT payload FROM curves WHERE id={id}"))).unwrap()
}

impl Spec {
    fn validate(&self) {
        assert!(self.p < (1u64 << 62) && is_prime_u64(self.p) && is_prime_u64(self.r));
        assert!(
            self.a < self.p
                && self.b < self.p
                && self.public_check_scalar > 0
                && self.public_check_scalar < self.r
        );
        assert!(
            self.runs > 0
                && self.trials > 0
                && self.trials <= i64::MAX as u64
                && self.runs <= i64::MAX as u64
                && !self.primes.is_empty()
        );
        self.seed_start.checked_add(self.runs).unwrap();
        assert!(self.primes.windows(2).all(|w| w[0] < w[1]));
        assert!(self
            .primes
            .iter()
            .all(|&l| l >= 3 && l < self.p && l != self.r && is_prime_u64(l)));
        let t = self.p as i128 + 1 - self.r as i128;
        let d = t * t - 4 * self.p as i128;
        assert!(
            d < 0 && d.rem_euclid(4) == 1,
            "ordinary maximal-order fixture required"
        );
        let factors = factor_u64((-d).try_into().unwrap());
        assert!(
            factors.iter().all(|&(p, e)| is_prime_u64(p) && e == 1),
            "fundamental Frobenius discriminant required"
        );
        assert_eq!(
            factors
                .iter()
                .fold(1u128, |n, &(p, e)| n * (p as u128).pow(e)),
            (-d) as u128,
            "factorization must reconstruct the exact discriminant"
        );
        let f = Field::new(&self.p.into()).unwrap();
        let n = make_node(self, &f, self.a, self.b, self.g, self.q, 0, None);
        check_points(self, &n);
        // A nonzero r-torsion point and this unique Hasse multiple prove #E=r.
        assert!(2 * (self.r as u128) > self.p as u128 + 1 + 2 * integer_sqrt(self.p as u128));
    }
    #[cfg(test)]
    fn start(&self) -> StartCurve {
        StartCurve {
            name: "disk-walk-source".into(),
            tag: "fp".into(),
            p: self.p.into(),
            a: self.a.into(),
            b: self.b.into(),
            order: self.r.into(),
            subgroup_order: self.r.into(),
            cofactor: 1u8.into(),
            gx: self.g[0].into(),
            gy: self.g[1].into(),
        }
    }
}
fn integer_sqrt(n: u128) -> u128 {
    let (mut lo, mut hi) = (0, n + 1);
    while hi - lo > 1 {
        let m = (lo + hi) / 2;
        if m * m <= n {
            lo = m
        } else {
            hi = m
        }
    }
    lo
}
fn make_node(
    s: &Spec,
    f: &Field,
    a: u64,
    b: u64,
    g: [u64; 2],
    q: [u64; 2],
    depth: u64,
    parent_edge: Option<i64>,
) -> Node {
    let m = Model {
        a: f.from_u64(a),
        b: f.from_u64(b),
    };
    let j = word(f, &m.j(f).expect("nonsingular model"));
    let (fr, cr) = ec1_records(
        &s.p.into(),
        &a.into(),
        &b.into(),
        &s.r.into(),
        &1u8.into(),
        (&g[0].into(), &g[1].into()),
    );
    let (ec1, uid) = ec1_identity(&s.p.into(), &fr, &cr, "fp");
    Node {
        j,
        a,
        b,
        g,
        q,
        depth,
        parent_edge,
        ec1,
        uid,
    }
}
fn fast(s: &Spec, n: &Node) -> FastCurve {
    FastCurve::new("", s.p, n.a, n.b, s.r)
}
fn check_points(s: &Spec, n: &Node) {
    let c = fast(s, n);
    let g = c.point(n.g[0], n.g[1]);
    let q = c.point(n.q[0], n.q[1]);
    assert!(n.g.iter().chain(&n.q).all(|&x| x < s.p));
    assert!(c.is_on_curve(g) && c.is_on_curve(q));
    assert!(c.mul(Some(g), s.r).is_none() && c.mul(Some(q), s.r).is_none());
    assert_eq!(c.mul(Some(g), s.public_check_scalar), Some(q));
}
fn open(s: &Spec, path: &str) -> Db {
    s.validate();
    let db = Db::open(path);
    db.exec("CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS curves(id INTEGER PRIMARY KEY,j INTEGER UNIQUE NOT NULL,payload TEXT NOT NULL,expanded INTEGER NOT NULL DEFAULT 0 CHECK(expanded IN(0,1)),measured INTEGER NOT NULL DEFAULT 0 CHECK(measured IN(0,1)));
        CREATE INDEX IF NOT EXISTS frontier ON curves(expanded,id);
        CREATE INDEX IF NOT EXISTS measurement_queue ON curves(measured,id);
        CREATE TABLE IF NOT EXISTS edges(id INTEGER PRIMARY KEY,source INTEGER NOT NULL REFERENCES curves(id),target INTEGER NOT NULL REFERENCES curves(id),payload TEXT NOT NULL);
        CREATE INDEX IF NOT EXISTS edge_sources ON edges(source);
        CREATE TABLE IF NOT EXISTS failures(id INTEGER PRIMARY KEY,source INTEGER NOT NULL REFERENCES curves(id),payload TEXT NOT NULL);
        CREATE INDEX IF NOT EXISTS failure_sources ON failures(source);
        CREATE TABLE IF NOT EXISTS runs(curve INTEGER NOT NULL REFERENCES curves(id),run INTEGER NOT NULL,payload TEXT NOT NULL,targets INTEGER NOT NULL,relations INTEGER NOT NULL,verified INTEGER NOT NULL,PRIMARY KEY(curve,run));");
    db.exec("BEGIN IMMEDIATE");
    let spec = encode(s);
    if db.int("SELECT count(*) FROM meta") == 0 {
        let f = Field::new(&s.p.into()).unwrap();
        let root = make_node(s, &f, s.a, s.b, s.g, s.q, 0, None);
        db.exec(&format!("INSERT INTO meta VALUES('schema','isogeny-disk-walk/v1'),('spec',{}); INSERT INTO curves(id,j,payload) VALUES(1,{},{});",quote(&spec),root.j,quote(&encode(&root))));
    } else {
        assert_eq!(
            db.text("SELECT value FROM meta WHERE key='schema'"),
            "isogeny-disk-walk/v1"
        );
        assert_eq!(
            db.text("SELECT value FROM meta WHERE key='spec'"),
            spec,
            "refusing a changed field, target, degree set or run policy"
        );
    }
    db.exec("COMMIT");
    db
}
struct Engine<'a> {
    s: &'a Spec,
    f: Field,
    phis: Vec<ModPoly>,
}
impl<'a> Engine<'a> {
    fn new(s: &'a Spec) -> Self {
        let f = Field::new(&s.p.into()).unwrap();
        let phis = s
            .primes
            .iter()
            .map(|&l| modpoly::modular_polynomial(&f, l).expect("native modular polynomial"))
            .collect();
        Self { s, f, phis }
    }
    fn advance(&self, db: &Db, target: i64) {
        let mut completed = db.int("SELECT count(*) FROM curves WHERE expanded=1");
        while completed < target {
            db.exec("BEGIN IMMEDIATE");
            let id = db.int("SELECT coalesce(min(id),0) FROM curves WHERE expanded=0");
            if id == 0 {
                db.exec("COMMIT");
                break;
            }
            let src = read_node(db, id);
            let m = model(&self.f, &src);
            for phi in &self.phis {
                for j2 in poly::roots(
                    &self.f,
                    &phi.at_x(&self.f, &self.f.from_u64(src.j)),
                    1 ^ phi.ell,
                ) {
                    let iso = match kernel::elkies_isogeny(&self.f, &m, phi, &j2) {
                        Ok(i) => i,
                        Err(e) => {
                            self.failure(db, id, phi.ell, word(&self.f, &j2), &format!("{e:?}"));
                            continue;
                        }
                    };
                    let j = word(&self.f, &j2);
                    let dst_id = db.int(&format!(
                        "SELECT coalesce((SELECT id FROM curves WHERE j={j}),0)"
                    ));
                    let target_model = if dst_id != 0 {
                        model(&self.f, &read_node(db, dst_id))
                    } else {
                        let Some((m, _)) = curve::canonical_model(&self.f, &iso.codomain) else {
                            self.failure(db, id, phi.ell, j, "SpecialJ");
                            continue;
                        };
                        m
                    };
                    let map =
                        KernelMap::new(&self.f, &m, &target_model, phi.ell, iso.kernel.clone());
                    let image_g = pair(map.eval(tuple(src.g)));
                    let image_q = pair(map.eval(tuple(src.q)));
                    let edge_id = db.int("SELECT coalesce(max(id),0)+1 FROM edges");
                    let dst = make_node(
                        self.s,
                        &self.f,
                        word(&self.f, &target_model.a),
                        word(&self.f, &target_model.b),
                        image_g,
                        image_q,
                        src.depth + 1,
                        Some(edge_id),
                    );
                    self.check_map(&src, &dst, &map);
                    let target = if dst_id == 0 {
                        db.exec(&format!(
                            "INSERT INTO curves(j,payload) VALUES({j},{})",
                            quote(&encode(&dst))
                        ));
                        db.int("SELECT last_insert_rowid()")
                    } else {
                        dst_id
                    };
                    let edge = Edge {
                        source: id,
                        target,
                        degree: phi.ell,
                        kernel: iso.kernel.iter().map(|x| word(&self.f, x)).collect(),
                        u2: word(&self.f, &map.u2),
                        u3: word(&self.f, &map.u3),
                        image_g,
                        image_q,
                    };
                    db.exec(&format!("INSERT INTO edges(id,source,target,payload) VALUES({edge_id},{id},{target},{})",quote(&encode(&edge))));
                }
            }
            db.exec(&format!("UPDATE curves SET expanded=1 WHERE id={id}"));
            #[cfg(test)]
            if std::env::var_os("ISOGENY_DISK_TEST_ABORT").is_some() {
                std::process::exit(86)
            }
            db.exec("COMMIT");
            completed += 1;
            println!("{}", json!({"kind":"expanded","curve":id}));
            io::stdout().flush().unwrap();
        }
    }
    fn failure(&self, db: &Db, id: i64, ell: u64, j: u64, error: &str) {
        let row = json!({"source":id,"degree":ell,"target_j":j.to_string(),"error":error});
        db.exec(&format!(
            "INSERT INTO failures(source,payload) VALUES({id},{})",
            quote(&encode(&row))
        ));
    }
    fn check_map(&self, src: &Node, dst: &Node, map: &KernelMap<'_>) {
        check_points(self.s, src);
        check_points(self.s, dst);
        let c = fast(self.s, src);
        let d = fast(self.s, dst);
        for k in [0, 1, 2, 3, 7, 17, 31, 63] {
            let a = c
                .mul(Some(c.point(src.g[0], src.g[1])), k)
                .map(|p| c.canonical(p));
            let b = d
                .mul(Some(d.point(dst.g[0], dst.g[1])), k)
                .map(|p| d.canonical(p));
            assert_eq!(map.eval(a), b, "scalar transport");
        }
        assert_eq!(map.eval(tuple(src.q)), tuple(dst.q));
    }
    fn verify(&self, db: &Db, replay_first: i64) -> Value {
        assert_eq!(db.text("PRAGMA integrity_check"), "ok");
        let mut foreign = db.prepare("PRAGMA foreign_key_check");
        assert!(!foreign.step());
        drop(foreign);
        let mut nodes = db.prepare("SELECT id,j,payload FROM curves ORDER BY id");
        while nodes.step() {
            let id = nodes.int(0);
            let n: Node = serde_json::from_str(&nodes.text(2)).unwrap();
            check_points(self.s, &n);
            assert_eq!(n.j, nodes.int(1) as u64);
            let rebuilt = make_node(self.s, &self.f, n.a, n.b, n.g, n.q, n.depth, n.parent_edge);
            assert_eq!(rebuilt, n);
            if id == 1 {
                assert_eq!(
                    (n.a, n.b, n.g, n.q),
                    (self.s.a, self.s.b, self.s.g, self.s.q)
                );
                assert_eq!(n.parent_edge, None)
            } else {
                let eid = n.parent_edge.expect("parent edge");
                let e: Edge = serde_json::from_str(
                    &db.text(&format!("SELECT payload FROM edges WHERE id={eid}")),
                )
                .unwrap();
                assert_eq!(e.target, id);
                assert_eq!(e.image_g, n.g);
                assert_eq!(e.image_q, n.q);
                assert_eq!(n.depth, read_node(db, e.source).depth + 1)
            }
        }
        let mut edges = db.prepare("SELECT source,target,payload FROM edges ORDER BY id");
        let mut count = 0;
        while edges.step() {
            let e: Edge = serde_json::from_str(&edges.text(2)).unwrap();
            assert_eq!((e.source, e.target), (edges.int(0), edges.int(1)));
            assert!(self.s.primes.contains(&e.degree));
            let src = read_node(db, e.source);
            let mut dst = read_node(db, e.target);
            let h = e.kernel.iter().map(|&x| self.f.from_u64(x)).collect();
            let map = KernelMap::new(
                &self.f,
                &model(&self.f, &src),
                &model(&self.f, &dst),
                e.degree,
                h,
            );
            assert_eq!(
                [word(&self.f, &map.u2), word(&self.f, &map.u3)],
                [e.u2, e.u3]
            );
            dst.g = e.image_g;
            dst.q = e.image_q;
            self.check_map(&src, &dst, &map);
            assert_eq!(map.eval(tuple(src.g)), tuple(dst.g));
            count += 1;
        }
        let mut sources = db.prepare("SELECT id FROM curves WHERE expanded=1 ORDER BY id");
        let mut expansions = 0;
        while sources.step() {
            let id = sources.int(0);
            let n = read_node(db, id);
            let mut expected = Vec::new();
            for phi in &self.phis {
                for j in poly::roots(
                    &self.f,
                    &phi.at_x(&self.f, &self.f.from_u64(n.j)),
                    1 ^ phi.ell,
                ) {
                    expected.push((phi.ell, word(&self.f, &j)));
                }
            }
            let mut got = Vec::new();
            let mut edge_rows = db.prepare(&format!("SELECT payload FROM edges WHERE source={id}"));
            while edge_rows.step() {
                let e: Edge = serde_json::from_str(&edge_rows.text(0)).unwrap();
                got.push((e.degree, read_node(db, e.target).j));
            }
            let mut failed = db.prepare(&format!("SELECT payload FROM failures WHERE source={id}"));
            while failed.step() {
                let e: Value = serde_json::from_str(&failed.text(0)).unwrap();
                got.push((
                    e["degree"].as_u64().unwrap(),
                    e["target_j"].as_str().unwrap().parse::<u64>().unwrap(),
                ));
            }
            expected.sort();
            got.sort();
            assert_eq!(expected, got, "expanded-source rational root roster");
            expansions += 1;
        }
        assert_eq!(
            db.int(
                "SELECT count(*) FROM edges e JOIN curves c ON e.source=c.id WHERE c.expanded=0"
            ),
            0
        );
        let mut rows = db.prepare(
            "SELECT curve,run,payload,targets,relations,verified FROM runs ORDER BY curve,run",
        );
        let mut replayed = 0;
        let mut relation_checks = 0;
        while rows.step() {
            let id = rows.int(0);
            let run = rows.int(1);
            let got: Run = serde_json::from_str(&rows.text(2)).unwrap();
            let n = read_node(db, id);
            assert!(run >= 0 && (run as u64) < self.s.runs);
            assert_eq!(got.seed, self.s.seed_start + run as u64);
            assert_eq!(
                (
                    got.targets as i64,
                    got.relations.len() as i64,
                    i64::from(got.scalar.is_some())
                ),
                (rows.int(3), rows.int(4), rows.int(5))
            );
            validate_run(self.s, &n, &got);
            relation_checks += got.relations.len();
            if id <= replay_first {
                assert_eq!(got, collect(self.s, &n, run as u64));
                replayed += 1
            }
        }
        assert_eq!(db.int(&format!("SELECT count(*) FROM curves c WHERE measured != ((SELECT count(*) FROM runs r WHERE r.curve=c.id)={})",self.s.runs)),0);
        assert_eq!(db.int("SELECT count(*) FROM curves c WHERE (SELECT count(*) FROM runs r WHERE r.curve=c.id)!=(SELECT coalesce(max(run)+1,0) FROM runs r WHERE r.curve=c.id)"),0);
        json!({"schema":"isogeny-disk-verify/v1","result":"PASS","curves":db.int("SELECT count(*) FROM curves"),"edges_reverified":count,"expanded_source_rosters_reverified":expansions,"scalar_transport_checks":count*8,"relations_readded":relation_checks,"collection_runs_replayed":replayed,"independent_host":false})
    }
}
fn collect(s: &Spec, n: &Node, run: u64) -> Run {
    check_points(s, n);
    let c = fast(s, n);
    let base = FastFactorBase::build(&c, 3);
    for &p in &base.pts {
        assert!(c.is_on_curve(p) && c.mul(Some(p), s.r).is_none())
    }
    let (rels, stats) = collect_relations(
        &c,
        c.point(n.g[0], n.g[1]),
        c.point(n.q[0], n.q[1]),
        &base,
        None,
        20,
        s.trials,
        s.seed_start + run,
    );
    let scalar = solve_relations(&c, 3, &rels);
    let row = Run {
        seed: s.seed_start + run,
        targets: stats.targets,
        relations: rels
            .into_iter()
            .map(|r| Relation {
                a: r.a,
                b: r.b,
                entries: r.entries,
            })
            .collect(),
        rejected: stats.rejected,
        scalar,
        status: if scalar.is_some() {
            "VERIFIED_SOLVE"
        } else if stats.targets == s.trials {
            "EXHAUSTED_TRIAL_CAP"
        } else {
            "UNRESOLVED_RANK"
        }
        .into(),
    };
    validate_run(s, n, &row);
    row
}
fn validate_run(s: &Spec, n: &Node, row: &Run) {
    assert_eq!(row.rejected, 0);
    assert!(row.targets <= s.trials);
    let c = fast(s, n);
    let base = FastFactorBase::build(&c, 3);
    let g = c.point(n.g[0], n.g[1]);
    let q = c.point(n.q[0], n.q[1]);
    let mut rels = Vec::new();
    for rel in &row.relations {
        let mut rhs = None;
        for &(i, sign) in &rel.entries {
            let p = base.pts[i as usize];
            let p = if sign < 0 { c.neg(p) } else { p };
            rhs = c.add(rhs, c.mul(Some(p), sign.unsigned_abs() as u64));
        }
        assert_eq!(c.combine(g, q, rel.a, rel.b), rhs);
        rels.push(FastRelation {
            a: rel.a,
            b: rel.b,
            entries: rel.entries.clone(),
        });
    }
    assert_eq!(solve_relations(&c, 3, &rels), row.scalar);
    if let Some(k) = row.scalar {
        assert_eq!(k, s.public_check_scalar);
        assert_eq!(c.mul(Some(g), k), Some(q));
        assert_eq!(row.status, "VERIFIED_SOLVE")
    } else if row.targets == s.trials {
        assert_eq!(row.status, "EXHAUSTED_TRIAL_CAP")
    } else {
        assert_eq!(row.status, "UNRESOLVED_RANK")
    }
}
fn measure(s: &Spec, db: &Db, first: i64) {
    loop {
        db.exec("BEGIN IMMEDIATE");
        let id = db.int(&format!(
            "SELECT coalesce(min(id),0) FROM curves WHERE measured=0 AND id<={first}"
        ));
        if id == 0 {
            db.exec("COMMIT");
            break;
        }
        let n = read_node(db, id);
        let run = db.int(&format!("SELECT count(*) FROM runs WHERE curve={id}"));
        assert!((run as u64) <= s.runs);
        if (run as u64) < s.runs {
            let row = collect(s, &n, run as u64);
            db.exec(&format!(
                "INSERT INTO runs VALUES({id},{run},{},{},{},{})",
                quote(&encode(&row)),
                row.targets,
                row.relations.len(),
                i32::from(row.scalar.is_some())
            ));
        }
        if run as u64 + 1 >= s.runs {
            db.exec(&format!("UPDATE curves SET measured=1 WHERE id={id}"))
        }
        db.exec("COMMIT");
    }
}
fn status(s: &Spec, db: &Db, path: &str) -> Value {
    let pending = db.int("SELECT count(*) FROM curves WHERE expanded=0");
    let failures = db.int("SELECT count(*) FROM failures");
    let rss = fs::read_to_string("/proc/self/status").ok().and_then(|x| {
        x.lines()
            .find(|l| l.starts_with("VmHWM:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|n| n.parse::<u64>().ok())
    });
    json!({"schema":"isogeny-disk-status/v1","discovered":db.int("SELECT count(*) FROM curves"),"expanded":db.int("SELECT count(*) FROM curves WHERE expanded=1"),"pending":pending,"edges":db.int("SELECT count(*) FROM edges"),"failures":failures,"measured_curves":db.int("SELECT count(*) FROM curves WHERE measured=1"),"runs":db.int("SELECT count(*) FROM runs"),"target_trials":db.int("SELECT coalesce(sum(targets),0) FROM runs"),"relations":db.int("SELECT coalesce(sum(relations),0) FROM runs"),"verified_solves":db.int("SELECT coalesce(sum(verified),0) FROM runs"),"runs_per_curve":s.runs,"sqlite_cache_bytes":db.cache_bytes(),"process_peak_rss_kib":rss,"database_bytes":fs::metadata(path).unwrap().len(),"selected_prime_component_closed":pending==0&&failures==0,"entire_horizontal_class_certified":false,"full_cost":null,"speedup":null})
}
fn dump(db: &Db, mut out: impl Write) {
    for (kind, query) in [
        ("meta", "SELECT key,value FROM meta ORDER BY key"),
        (
            "curve",
            "SELECT id,payload,expanded,measured FROM curves ORDER BY id",
        ),
        ("edge", "SELECT id,payload FROM edges ORDER BY id"),
        ("failure", "SELECT id,payload FROM failures ORDER BY id"),
        (
            "run",
            "SELECT curve,run,payload FROM runs ORDER BY curve,run",
        ),
    ] {
        let mut q = db.prepare(query);
        while q.step() {
            let v = if kind == "meta" {
                json!({"kind":kind,"key":q.text(0),"value":q.text(1)})
            } else if kind == "run" {
                json!({"kind":kind,"curve":q.int(0),"run":q.int(1),"payload":serde_json::from_str::<Value>(&q.text(2)).unwrap()})
            } else if kind == "curve" {
                json!({"kind":kind,"id":q.int(0),"payload":serde_json::from_str::<Value>(&q.text(1)).unwrap(),"expanded":q.int(2),"measured":q.int(3)})
            } else {
                json!({"kind":kind,"id":q.int(0),"payload":serde_json::from_str::<Value>(&q.text(1)).unwrap()})
            };
            writeln!(out, "{v}").unwrap();
        }
    }
}
pub fn main() {
    let a: Vec<String> = std::env::args().collect();
    assert!(a.len()>=4,"usage: isogeny_disk_walk init|advance|measure|verify|status|dump|snapshot SPEC.json DB [cumulative-limit|replay-first-curves|snapshot-path]");
    let s: Spec = serde_json::from_str(&fs::read_to_string(&a[2]).unwrap()).unwrap();
    let db = open(&s, &a[3]);
    let limit = || {
        a.get(4)
            .expect("limit")
            .parse::<i64>()
            .expect("nonnegative integer limit")
    };
    match a[1].as_str() {
        "init" => (),
        "advance" => {
            assert!(limit() >= 0);
            Engine::new(&s).advance(&db, limit())
        }
        "measure" => {
            assert!(limit() >= 0);
            measure(&s, &db, limit())
        }
        "verify" => {
            let e = Engine::new(&s);
            println!(
                "{}",
                e.verify(&db, a.get(4).map(|v| v.parse().unwrap()).unwrap_or(0))
            )
        }
        "snapshot" => {
            let dest = a.get(4).expect("new snapshot path");
            assert!(!Path::new(dest).exists());
            db.exec(&format!("VACUUM INTO {}", quote(dest)))
        }
        "dump" => {
            dump(&db, io::BufWriter::new(io::stdout()));
            return;
        }
        "status" => (),
        _ => panic!("unknown command"),
    }
    println!("{}", status(&s, &db, &a[3]));
}

#[cfg(test)]
#[path = "isogeny_disk/tests.rs"]
mod tests;

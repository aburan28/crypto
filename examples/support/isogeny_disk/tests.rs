use super::*;
use crypto_lib::cryptanalysis::isogeny_walk::walk::{Walk, WalkConfig};
use std::collections::BTreeSet;
use std::sync::atomic::{AtomicU64, Ordering};
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Temp(std::path::PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "isogeny-disk-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::SeqCst)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
    fn file(&self, n: &str) -> String {
        self.0.join(n).to_str().unwrap().into()
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn count(p: u64, a: u64, b: u64) -> u64 {
    let mut squares = vec![false; p as usize];
    for y in 0..p {
        squares[(y * y % p) as usize] = true
    }
    1 + (0..p)
        .map(|x| {
            let v = (x * x % p * x + a * x + b) % p;
            if v == 0 {
                1
            } else if squares[v as usize] {
                2
            } else {
                0
            }
        })
        .sum::<u64>()
}
fn fixture(p: u64, r: u64, ell: u64) -> Spec {
    for a in 1..p {
        for b in 1..p {
            if count(p, a, b) != r {
                continue;
            }
            let c = FastCurve::new("", p, a, b, r);
            for x in 0..p {
                for y in 1..p {
                    let g = c.point(x, y);
                    if !c.is_on_curve(g) {
                        continue;
                    }
                    let (qx, qy) = c.canonical(c.mul(Some(g), 7).unwrap());
                    return Spec {
                        p,
                        r,
                        a,
                        b,
                        g: [x, y],
                        q: [qx, qy],
                        public_check_scalar: 7,
                        primes: vec![ell],
                        runs: 2,
                        seed_start: 2026101600,
                        trials: 100000,
                    };
                }
            }
        }
    }
    panic!("no fixture")
}
fn all_j(s: &Spec) -> BTreeSet<u64> {
    let f = Field::new(&s.p.into()).unwrap();
    let mut out = BTreeSet::new();
    for a in 1..s.p {
        for b in 1..s.p {
            if count(s.p, a, b) == s.r {
                let m = Model {
                    a: f.from_u64(a),
                    b: f.from_u64(b),
                };
                out.insert(word(&f, &m.j(&f).unwrap()));
            }
        }
    }
    out
}
fn graph(db: &Db) -> BTreeSet<(u64, u64, Vec<u64>, u64)> {
    let mut set = BTreeSet::new();
    let mut q = db.prepare("SELECT payload FROM edges ORDER BY id");
    while q.step() {
        let e: Edge = serde_json::from_str(&q.text(0)).unwrap();
        set.insert((
            read_node(db, e.source).j,
            e.degree,
            e.kernel,
            read_node(db, e.target).j,
        ));
    }
    set
}
#[test]
fn complete_small_classes_match_independent_point_counts_and_existing_walker() {
    for (p, r, l) in [(101, 97, 5), (127, 113, 7)] {
        let t = Temp::new();
        let s = fixture(p, r, l);
        let db = open(&s, &t.file("walk.db"));
        let e = Engine::new(&s);
        e.advance(&db, 100);
        let mut w = Walk::new(
            s.start(),
            WalkConfig {
                primes: s.primes.clone(),
                max_curves: 256,
                max_depth: usize::MAX,
                audit_points: 2,
                ..Default::default()
            },
        )
        .unwrap();
        w.run();
        let wanted: BTreeSet<_> = w
            .edges
            .iter()
            .map(|x| {
                (
                    word(&w.field, &w.nodes[x.source].j),
                    x.ell,
                    x.kernel
                        .iter()
                        .map(|x| word(&w.field, x))
                        .collect::<Vec<_>>(),
                    word(&w.field, &w.nodes[x.target].j),
                )
            })
            .collect();
        assert_eq!(graph(&db), wanted);
        assert!(w.failures.is_empty());
        assert_eq!(db.int("SELECT count(*) FROM failures"), 0);
        let mut js = BTreeSet::new();
        let mut rows = db.prepare("SELECT j FROM curves");
        while rows.step() {
            js.insert(rows.int(0) as u64);
        }
        assert_eq!(js, all_j(&s));
        assert_eq!(js.len(), 3);
        measure(&s, &db, 100);
        assert_eq!(db.int("SELECT count(*) FROM runs"), 6);
        assert_eq!(db.int("SELECT sum(verified) FROM runs"), 6);
        let receipt = e.verify(&db, 100);
        assert_eq!(receipt["result"], "PASS");
        assert_eq!(receipt["collection_runs_replayed"], 6);
        assert!(receipt["relations_readded"].as_u64().unwrap() > 0);
    }
}
#[test]
fn crash_child() {
    if let Ok(dir) = std::env::var("ISOGENY_DISK_CRASH_CHILD") {
        let s: Spec =
            serde_json::from_str(&fs::read_to_string(format!("{dir}/spec.json")).unwrap()).unwrap();
        let db = open(&s, &format!("{dir}/resume.db"));
        Engine::new(&s).advance(&db, 100);
        panic!("did not interrupt before commit")
    }
}
#[test]
fn abrupt_exit_rolls_back_and_resume_is_identical_to_uninterrupted() {
    let t = Temp::new();
    let s = fixture(101, 97, 5);
    fs::write(t.file("spec.json"), encode(&s)).unwrap();
    drop(open(&s, &t.file("resume.db")));
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "isogeny_disk::tests::crash_child", "--nocapture"])
        .env("ISOGENY_DISK_CRASH_CHILD", &t.0)
        .env("ISOGENY_DISK_TEST_ABORT", "1")
        .status()
        .unwrap();
    assert_eq!(result.code(), Some(86));
    let db = open(&s, &t.file("resume.db"));
    assert_eq!(db.int("SELECT count(*) FROM curves"), 1);
    assert_eq!(db.int("SELECT count(*) FROM edges"), 0);
    assert_eq!(db.int("SELECT sum(expanded) FROM curves"), 0);
    let e = Engine::new(&s);
    e.advance(&db, 1);
    e.advance(&db, 100);
    measure(&s, &db, 1);
    measure(&s, &db, 100);
    let mut a = Vec::new();
    dump(&db, &mut a);
    let other = open(&s, &t.file("direct.db"));
    e.advance(&other, 100);
    measure(&s, &other, 100);
    let mut b = Vec::new();
    dump(&other, &mut b);
    assert_eq!(a, b);
    let snap = t.file("snapshot.db");
    db.exec(&format!("VACUUM INTO {}", quote(&snap)));
    let restored = open(&s, &snap);
    assert_eq!(e.verify(&restored, 100)["result"], "PASS");
}
#[test]
fn wrong_policy_and_corrupt_j_index_fail_closed() {
    let t = Temp::new();
    let s = fixture(101, 97, 5);
    drop(open(&s, &t.file("db")));
    let mut changed = s.clone();
    changed.trials += 1;
    assert!(std::panic::catch_unwind(|| open(&changed, &t.file("db"))).is_err());
    let db = open(&s, &t.file("db"));
    db.exec("UPDATE curves SET j=j+1 WHERE id=1");
    let e = Engine::new(&s);
    assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| e.verify(&db, 0))).is_err());
}

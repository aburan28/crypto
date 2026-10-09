//! Installed isogeny construction, screening, and independent certificate replay.
#![allow(dead_code)]
mod constructor;
#[path = "../../../src/cryptanalysis/isogeny_walk/curve.rs"]
mod curve;
#[path = "../../../src/cryptanalysis/isogeny_walk/field.rs"]
mod field;
#[path = "../../../src/hash/sha256.rs"]
mod hash_sha256;
#[path = "../../../research/p192_p224_large_degree_isogenies_20261009/verification_kernel.rs"]
mod kernel;
#[path = "../../../research/p192_p224_large_degree_isogenies_20261009/map_identity.rs"]
mod map_identity;
#[path = "../../../src/cryptanalysis/isogeny_walk/modpoly.rs"]
mod modpoly;
#[cfg(test)]
#[path = "../../../src/cryptanalysis/modular_polynomial.rs"]
pub mod modular_polynomial;
#[path = "../../../research/p192_p224_large_degree_isogenies_20261009/verification_poly.rs"]
mod poly;
mod replay;
mod supervisor;
#[cfg(test)]
pub mod cryptanalysis {
    pub use crate::modular_polynomial;
}

use isogeny_algos::{
    api::{self, PrimeCurve},
    field::{is_prime, Zp},
    icv1,
    int::Int,
    json::Json,
};
use serde_json::{json, Value};
use std::{
    collections::HashMap,
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

const REGISTRY: &str = include_str!("../../../docs/curves/registry.json");
const KERNEL_PROTOCOL: &[u8] = include_bytes!(
    "../../../research/p192_p224_large_degree_isogenies_20261009/KERNEL_PROTOCOL.md"
);
const GENERAL_PROTOCOL: &[u8] =
    include_bytes!("../../../research/p192_p224_large_degree_isogenies_20261009/PROTOCOL.md");
const HELP: &str = "isogeny — public-curve isogeny search and independent verification

  isogeny search --curve p224 --ell 1471 --out RESULTS [--timeout 1800]
  isogeny search --curve p192 --ell 10453 --out RESULTS [--timeout 1800]
  isogeny verify --input RESULTS
  isogeny screen --curve p192 --from 1010 --to 65537 [--max-order 8]
  isogeny version

search options:
  --method auto|kernel|modpoly  auto selects the validated kernel route for
                              p224/1471 and p192/10453; otherwise full enumeration
  --construct-only            record construction without independent replay
  --timeout SECONDS           construction time limit (default 1800; RSS cap 8 GiB)
  --out DIRECTORY             fresh output directory; existing paths are refused

Search normally verifies the result independently before reporting success.
The kernel route covers one Frobenius eigenline of two (coverage PARTIAL).
The modular-polynomial route enumerates both lines on a split prime.
The published source orders, seed 1, field arithmetic, and standard generators
are built into this executable; no checkout, Cargo, or companion program is needed.
Exit: 0 passed; 1 failed; 2 usage; 3 timeout/resource limit/incomplete construction.
";

fn sha256_hex(bytes: &[u8]) -> String {
    isogeny_algos::sha256::sha256_hex(bytes)
}
fn standard_order(name: &str) -> Int {
    supervisor::order(name)
}
fn provenance() -> Value {
    json!({"name":"isogeny", "version":env!("CARGO_PKG_VERSION"),
        "git_commit":env!("ISOGENY_GIT_COMMIT"), "git_dirty":env!("ISOGENY_GIT_DIRTY"),
        "standards_sha256":sha256_hex(REGISTRY.as_bytes())})
}

struct Options(HashMap<String, String>);
impl Options {
    fn parse(args: &[String], allowed: &[&str]) -> Result<Self, String> {
        let mut values = HashMap::new();
        let mut i = 0;
        while i < args.len() {
            let key = args[i]
                .strip_prefix("--")
                .ok_or_else(|| format!("unexpected argument {}", args[i]))?;
            if !allowed.contains(&key) {
                return Err(format!("unknown option --{key}"));
            }
            let value = if key == "construct-only" {
                "true".to_owned()
            } else {
                i += 1;
                args.get(i)
                    .filter(|s| !s.starts_with("--"))
                    .ok_or_else(|| format!("--{key} needs a value"))?
                    .clone()
            };
            if values.insert(key.to_string(), value).is_some() {
                return Err(format!("duplicate --{key}"));
            }
            i += 1;
        }
        Ok(Self(values))
    }
    fn get(&self, key: &str) -> Option<&str> {
        self.0.get(key).map(String::as_str)
    }
    fn need(&self, key: &str) -> Result<&str, String> {
        self.get(key).ok_or_else(|| format!("--{key} is required"))
    }
    fn number(&self, key: &str, default: Option<u64>) -> Result<u64, String> {
        match self.get(key) {
            Some(s) => s
                .parse()
                .map_err(|_| format!("--{key} must be a nonnegative integer")),
            None => default.ok_or_else(|| format!("--{key} is required")),
        }
    }
    fn curve(&self) -> Result<&str, String> {
        let curve = self.need("curve")?;
        if matches!(curve, "p192" | "p224") {
            Ok(curve)
        } else {
            Err("--curve must be p192 or p224".into())
        }
    }
}

fn eigenvalues(name: &str, ell: u64) -> Option<[u64; 2]> {
    let c = PrimeCurve::preset(name).unwrap();
    let t = &(&c.p + &Int::one()) - &standard_order(name);
    let tm = t.mod_u64(ell);
    let pm = c.p.mod_u64(ell);
    let disc = (tm * tm + ell - (4 * pm) % ell) % ell;
    if disc == 0 {
        return None;
    }
    let field = Zp::new(ell);
    let sqrt = field.sqrt(disc)?;
    Some([
        (tm + sqrt) % ell * ((ell + 1) / 2) % ell,
        (tm + ell - sqrt) % ell * ((ell + 1) / 2) % ell,
    ])
}
fn multiplicative_order(a: u64, p: u64) -> u64 {
    let mut n = p - 1;
    let mut residual = n;
    let mut d = 2;
    while d * d <= residual {
        if residual % d == 0 {
            while n % d == 0 && supervisor::pow(a, n / d, p) == 1 {
                n /= d;
            }
            while residual % d == 0 {
                residual /= d;
            }
        }
        d += 1;
    }
    if residual > 1 {
        while n % residual == 0 && supervisor::pow(a, n / residual, p) == 1 {
            n /= residual;
        }
    }
    n
}

fn general_worker(name: &str, ell: u64) -> Result<Json, String> {
    let c = PrimeCurve::preset(name).unwrap();
    let n = standard_order(name);
    let identity = |c: &PrimeCurve| {
        let id = icv1::prime(&c.p, &c.a, &c.b, &n).ok_or("curve identity could not be formed")?;
        Ok::<_, String>(Json::obj(vec![
            ("slug", Json::str(id.slug)),
            ("icv1", Json::str(id.icv1)),
            ("model_json", Json::str(id.model_json)),
        ]))
    };
    let records = api::isogenies(&c, ell, 1)?;
    if records.len() != 2 || records.iter().any(|r| !r.verified()) {
        return Err("complete two-eigenline construction did not pass".into());
    }
    let array = |v: &[Int]| Json::Arr(v.iter().map(Json::str).collect());
    let mut maps = vec![];
    for r in records {
        let target = PrimeCurve::new(c.p.clone(), r.codomain.0, r.codomain.1)?;
        maps.push(Json::obj(vec![
            ("verified", Json::Bool(true)),
            ("kernel", array(&r.kernel)),
            (
                "x_map",
                Json::obj(vec![("num", array(&r.map_num)), ("den", array(&r.map_den))]),
            ),
            (
                "codomain",
                Json::obj(vec![
                    ("a", Json::str(&target.a)),
                    ("b", Json::str(&target.b)),
                    ("j", Json::str(r.j_codomain)),
                    ("icv1", identity(&target)?),
                ]),
            ),
        ]));
    }
    Ok(Json::obj(vec![
        ("status", Json::str("PASS")),
        ("scope", Json::str("all_frobenius_eigenlines")),
        ("degree_coverage", Json::str("COMPLETE")),
        (
            "curve",
            Json::obj(vec![
                ("preset", Json::str(name)),
                ("p", Json::str(&c.p)),
                ("a", Json::str(&c.a)),
                ("b", Json::str(&c.b)),
            ]),
        ),
        (
            "result",
            Json::obj(vec![
                ("ell", Json::Num(ell as i64)),
                ("source_icv1", identity(&c)?),
                ("isogenies", Json::Arr(maps)),
            ]),
        ),
    ]))
}

fn execute(args: &[String]) -> Result<(Value, i32), String> {
    let command = args.first().map(String::as_str).unwrap_or("help");
    match command {
        "version" | "--version" | "-V" => {
            Ok((json!({"schema":"isogeny-tool/v1","tool":provenance()}), 0))
        }
        "screen" => {
            let o = Options::parse(&args[1..], &["curve", "from", "to", "max-order"])?;
            let name = o.curve()?;
            let lower = o.number("from", Some(1010))?;
            let upper = o.number("to", Some(65537))?;
            let max_order = o.number("max-order", Some(8))?;
            if lower < 3 || upper < lower || upper > 1_000_000 || max_order == 0 {
                return Err(
                    "screen requires 3 <= from <= to <= 1000000 and positive max-order".into(),
                );
            }
            let mut candidates = vec![];
            for ell in lower..=upper {
                if !is_prime(ell) {
                    continue;
                }
                if let Some(roots) = eigenvalues(name, ell) {
                    let orders = roots.map(|a| multiplicative_order(a, ell));
                    if *orders.iter().min().unwrap() <= max_order {
                        candidates.push(json!({"curve":name,"ell":ell,"eigenvalues":roots,"eigenvalue_orders":orders}));
                    }
                }
            }
            Ok((
                json!({"schema":"isogeny-screen/v1","tool":provenance(),"lower":lower,"upper":upper,"maximum_selected_order":max_order,"candidates":candidates,"status":"PASS"}),
                0,
            ))
        }
        "search" => {
            let o = Options::parse(
                &args[1..],
                &["curve", "ell", "out", "timeout", "method", "construct-only"],
            )?;
            let name = o.curve()?;
            let ell = o.number("ell", None)?;
            let timeout = o.number("timeout", Some(1800))?;
            if !(3..=1_000_000).contains(&ell) || !is_prime(ell) || timeout == 0 {
                return Err(
                    "search needs an odd prime degree <= 1000000 and a positive timeout".into(),
                );
            }
            let supported = matches!((name, ell), ("p192", 10453) | ("p224", 1471));
            let single = match o.get("method").unwrap_or("auto") {
                "auto" => supported,
                "kernel" if supported => true,
                "kernel" => {
                    return Err(
                        "kernel construction currently supports p192/10453 and p224/1471".into(),
                    )
                }
                "modpoly" => false,
                _ => return Err("--method must be auto, kernel or modpoly".into()),
            };
            if eigenvalues(name, ell).is_none() {
                return Err(
                    "degree is inert or has a repeated Frobenius eigenvalue; select a split prime"
                        .into(),
                );
            }
            let out = PathBuf::from(o.need("out")?);
            if out.exists() {
                return Err("--out must be a fresh directory".into());
            }
            if let Some(parent) = out.parent().filter(|p| !p.as_os_str().is_empty()) {
                fs::create_dir_all(parent).map_err(|e| e.to_string())?;
            }
            fs::create_dir(&out).map_err(|e| e.to_string())?;
            let dir = out.join(name);
            fs::create_dir(&dir).map_err(|e| e.to_string())?;
            let executable = std::env::current_exe().map_err(|e| e.to_string())?;
            let worker_args = [
                if single { "kernel-first" } else { "isogenies" }.into(),
                "--curve".into(),
                name.into(),
                "--ell".into(),
                ell.to_string(),
                "--order".into(),
                standard_order(name).to_string(),
                "--seed".into(),
                "1".into(),
            ];
            let receipt = supervisor::run(
                &executable,
                &dir,
                &worker_args,
                &format!("ell-{ell}"),
                timeout,
            );
            let receipt: Value = serde_json::from_str(&receipt.dump()).unwrap();
            let summary = json!({"schema":"large-degree-isogeny-search/v1","source_commit":env!("ISOGENY_GIT_COMMIT"),"tool":provenance(),
                "executable_sha256":sha256_hex(&fs::read(&executable).map_err(|e| e.to_string())?),
                "protocol_sha256":sha256_hex(if single {KERNEL_PROTOCOL} else {GENERAL_PROTOCOL}),"timeout_seconds":timeout,
                "attempts":[{"curve":name,"ell":ell,"receipt":receipt}]});
            fs::write(
                out.join("search.json"),
                format!("{}\n", serde_json::to_string_pretty(&summary).unwrap()),
            )
            .map_err(|e| e.to_string())?;
            if receipt["status"] != "PASS" {
                let status = receipt["status"].as_str().unwrap_or("FAIL");
                return Ok((
                    json!({"status":status,"output":out,"construction":receipt,"tool":provenance()}),
                    if status == "FAIL" { 1 } else { 3 },
                ));
            }
            if o.get("construct-only").is_some() {
                Ok((
                    json!({"status":"CONSTRUCTED","independent_replay":"NOT_RUN","degree_coverage":if single {"PARTIAL"} else {"COMPLETE"},"output":out,"tool":provenance()}),
                    0,
                ))
            } else {
                let start = Instant::now();
                let certificate = replay::verify(&out, true);
                Ok((
                    json!({"status":"PASS","output":out,"degree_coverage":certificate["degree_coverage"],"verification_elapsed_ms":start.elapsed().as_millis(),"certificate":certificate,"tool":provenance()}),
                    0,
                ))
            }
        }
        "verify" => {
            let o = Options::parse(&args[1..], &["input"])?;
            let certificate = replay::verify(Path::new(o.need("input")?), false);
            Ok((
                json!({"status":"PASS","certificate":certificate,"tool":provenance()}),
                0,
            ))
        }
        "kernel-first" | "isogenies" => {
            let o = Options::parse(&args[1..], &["curve", "ell", "order", "seed"])?;
            let name = o.curve()?;
            let ell = o.number("ell", None)?;
            if o.need("order")? != standard_order(name).to_string() || o.need("seed")? != "1" {
                return Err("worker requires the published order and seed 1".into());
            }
            let record = if command == "kernel-first" {
                constructor::run(name, ell)
            } else {
                general_worker(name, ell)?
            };
            // Preserve the native record's coefficient ordering and exact bytes.
            println!("{}", record.dump());
            std::process::exit(0)
        }
        _ => Err(format!("unknown command {command}; run isogeny --help")),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty()
        || args
            .iter()
            .any(|a| matches!(a.as_str(), "help" | "--help" | "-h"))
    {
        println!("{HELP}");
        return;
    }
    let result = std::panic::catch_unwind(|| execute(&args));
    let (record, code) = match result {
        Ok(Ok(result)) => result,
        Ok(Err(error)) => (json!({"status":"USAGE","error":error}), 2),
        Err(_) => (
            json!({"status":"FAIL","error":"construction or independent verification failed; inspect retained run records"}),
            1,
        ),
    };
    println!("{}", serde_json::to_string_pretty(&record).unwrap());
    std::process::exit(code);
}

#[cfg(test)]
mod cli_tests {
    use super::*;
    #[test]
    fn known_small_extension_candidates_have_exact_eigenvalue_orders() {
        for (name, ell, root, order) in [("p192", 10453, 270, 3), ("p224", 1471, 554, 5)] {
            let roots = eigenvalues(name, ell).unwrap();
            assert!(roots.contains(&root));
            assert_eq!(multiplicative_order(root, ell), order);
        }
    }
    #[test]
    fn usage_failures_precede_output_creation() {
        for args in [
            vec![
                "search",
                "--curve",
                "p256",
                "--ell",
                "1471",
                "--out",
                "never-created",
            ],
            vec![
                "search",
                "--curve",
                "p224",
                "--ell",
                "1471",
                "--method",
                "wrong",
                "--out",
                "never-created",
            ],
        ] {
            let args: Vec<_> = args.into_iter().map(str::to_owned).collect();
            assert!(execute(&args).is_err());
        }
    }
}

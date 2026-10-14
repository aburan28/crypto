//! Standalone, reproducible BMSS benchmark. Counting and timing are separate executions.
use isogeny_algos::curve::*;
use isogeny_algos::field::{Field, Rng};
use isogeny_algos::find::bmss::{self, Method};
use isogeny_algos::fpm::FpM;
use isogeny_algos::json::Json;
use isogeny_algos::kernel::kohel::kohel;
use isogeny_algos::sha256::sha256_hex;
use isogeny_algos::testdata::supersingular_workload;
use std::fs::OpenOptions;
use std::hint::black_box;
use std::io::Write;
use std::sync::{
    atomic::{AtomicU64, Ordering},
    Arc,
};
use std::time::Instant;

#[path = "perf_campaign/mod.rs"]
mod campaign;

#[derive(Clone)]
struct Counted<F> {
    inner: F,
    counts: Arc<[AtomicU64; 4]>,
}
impl<F: Field> Counted<F> {
    fn new(inner: F) -> Self {
        Self {
            inner,
            counts: Arc::new(std::array::from_fn(|_| AtomicU64::new(0))),
        }
    }
    fn inc(&self, i: usize) {
        self.counts[i].fetch_add(1, Ordering::Relaxed);
    }
    fn snapshot(&self) -> Json {
        Json::obj(
            ["add_sub", "mul", "inv", "neg"]
                .into_iter()
                .enumerate()
                .map(|(i, k)| (k, Json::str(self.counts[i].load(Ordering::Relaxed))))
                .collect(),
        )
    }
}
// Counts use the generic Field implementation, not architecture-specific convolution
// instructions. They are an abstract arithmetic model; timing uses the original field.
impl<F: Field> Field for Counted<F> {
    type E = F::E;
    fn zero(&self) -> Self::E {
        self.inner.zero()
    }
    fn one(&self) -> Self::E {
        self.inner.one()
    }
    fn add(&self, a: Self::E, b: Self::E) -> Self::E {
        self.inc(0);
        self.inner.add(a, b)
    }
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E {
        self.inc(0);
        self.inner.sub(a, b)
    }
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E {
        self.inc(1);
        self.inner.mul(a, b)
    }
    fn inv(&self, a: Self::E) -> Self::E {
        self.inc(2);
        self.inner.inv(a)
    }
    fn neg(&self, a: Self::E) -> Self::E {
        self.inc(3);
        self.inner.neg(a)
    }
    fn from_u64(&self, n: u64) -> Self::E {
        self.inner.from_u64(n)
    }
    fn char(&self) -> u64 {
        self.inner.char()
    }
    fn characteristic_exceeds(&self, n: u64) -> bool {
        self.inner.characteristic_exceeds(n)
    }
    fn size(&self) -> u128 {
        self.inner.size()
    }
    fn q(&self) -> isogeny_algos::bigint::Big {
        self.inner.q()
    }
    fn random(&self, rng: &mut Rng) -> Self::E {
        self.inner.random(rng)
    }
    fn karatsuba_threshold(&self) -> usize {
        self.inner.karatsuba_threshold()
    }
}

fn command(program: &str, args: &[&str]) -> String {
    std::process::Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_owned())
        .unwrap_or_else(|| "unavailable".into())
}
fn host() -> Json {
    let cpu = if cfg!(target_os = "macos") {
        command("sysctl", &["-n", "machdep.cpu.brand_string"])
    } else {
        std::fs::read_to_string("/proc/cpuinfo")
            .map(|s| {
                s.lines()
                    .filter(|line| {
                        [
                            "model name",
                            "vendor_id",
                            "CPU implementer",
                            "CPU part",
                            "Hardware",
                            "flags",
                            "Features",
                        ]
                        .iter()
                        .any(|key| line.starts_with(key))
                    })
                    .collect::<std::collections::BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>()
                    .join("\n")
            })
            .unwrap_or_else(|_| "unavailable".into())
    };
    Json::obj(vec![
        ("os", Json::str(std::env::consts::OS)),
        ("arch", Json::str(std::env::consts::ARCH)),
        ("cpu", Json::str(cpu)),
        (
            "logical_cpus",
            Json::str(
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(0),
            ),
        ),
        ("kernel", Json::str(command("uname", &["-a"]))),
        (
            "isolation",
            Json::str("unreserved; wall samples are diagnostic and excluded from speedup claims"),
        ),
    ])
}
fn elements<const N: usize>(f: &FpM<N>, xs: &[[u64; N]]) -> Json {
    Json::Arr(xs.iter().map(|x| Json::str(f.to_big(x).to_dec())).collect())
}
fn map_record<const N: usize>(f: &FpM<N>, iso: &RatIsogeny<impl Field<E = [u64; N]>>) -> Json {
    Json::obj(vec![
        ("degree", Json::str(iso.deg)),
        ("a", Json::str(f.to_big(&iso.cod.a).to_dec())),
        ("b", Json::str(f.to_big(&iso.cod.b).to_dec())),
        ("kernel", elements(f, &iso.ker)),
        ("numerator", elements(f, &iso.num)),
        ("denominator", elements(f, &iso.den)),
    ])
}
fn csv<T: std::str::FromStr>(s: &str) -> Result<Vec<T>, String> {
    let v: Vec<T> = s
        .split(',')
        .map(|x| x.parse().map_err(|_| format!("invalid list: {s}")))
        .collect::<Result<_, _>>()?;
    if v.is_empty() {
        return Err("empty list".into());
    }
    Ok(v)
}
fn run<const N: usize>(
    out: &mut std::fs::File,
    bits: usize,
    ells: &[u64],
    samples: usize,
    timing: bool,
) -> Result<(), String> {
    let mut primes = vec![3];
    primes.extend_from_slice(ells);
    primes.sort();
    primes.dedup();
    let mut rng = Rng::new(0x1506_2026 + bits as u64);
    let w = supersingular_workload::<N>(&primes, bits, &mut rng);
    for &ell in ells {
        let p = w.pts.iter().find(|(l, _)| *l == ell).unwrap().1;
        let (g, _) = kernel_poly_from_point_f(&w.f, &w.e, &p, ell);
        let reference = kohel(&w.f, &w.e, &g, ell);
        let expected = map_record(&w.f, &reference).dump();
        let sigma = w.f.neg(w.f.mul(w.f.from_u64(2), g[g.len() - 2]));
        let input = Json::obj(vec![
            ("p", Json::str(w.f.q().to_dec())),
            ("a", Json::str(w.f.to_big(&w.e.a).to_dec())),
            ("b", Json::str(w.f.to_big(&w.e.b).to_dec())),
            (
                "codomain_a",
                Json::str(w.f.to_big(&reference.cod.a).to_dec()),
            ),
            (
                "codomain_b",
                Json::str(w.f.to_big(&reference.cod.b).to_dec()),
            ),
            ("ell", Json::str(ell)),
            ("sigma", Json::str(w.f.to_big(&sigma).to_dec())),
            ("reference_kernel", elements(&w.f, &g)),
        ]);
        let input_sha = sha256_hex(input.dump().as_bytes());
        let checks: Vec<_> = (0..4)
            .map(|_| random_point_f(&w.f, &w.e, &mut rng))
            .collect();
        for method in Method::ALL {
            eprintln!("{bits}-bit ell={ell} {}", method.name());
            let s = method.needs_sigma().then_some(sigma);
            let counted = Counted::new(w.f.clone());
            let actual = bmss::isogeny(&counted, method, &w.e, &reference.cod, ell as usize, s)
                .ok_or_else(|| format!("{} returned no map", method.name()))?;
            let counts = counted.snapshot();
            let output = map_record(&w.f, &actual).dump();
            if output != expected {
                return Err(format!("{} changed the exact map", method.name()));
            }
            let verified = bmss::isogeny(&w.f, method, &w.e, &reference.cod, ell as usize, s)
                .ok_or_else(|| format!("{} returned no native map", method.name()))?;
            if map_record(&w.f, &verified).dump() != expected
                || checks
                    .iter()
                    .any(|q| verified.eval(&w.f, q) != reference.eval(&w.f, q))
            {
                return Err(format!(
                    "{} native output or point transport mismatch",
                    method.name()
                ));
            }
            let mut raw = vec![];
            if timing {
                black_box(bmss::isogeny(
                    &w.f,
                    method,
                    &w.e,
                    &reference.cod,
                    ell as usize,
                    s,
                ));
                for _ in 0..samples {
                    let start = Instant::now();
                    let iso = black_box(bmss::isogeny(
                        black_box(&w.f),
                        method,
                        black_box(&w.e),
                        black_box(&reference.cod),
                        ell as usize,
                        s,
                    ));
                    let ns = start.elapsed().as_nanos();
                    let iso = iso.ok_or("timed call returned no map")?;
                    if map_record(&w.f, &iso).dump() != expected {
                        return Err("timed output mismatch".into());
                    }
                    raw.push(Json::str(ns));
                }
            }
            let row = Json::obj(vec![
                ("schema", Json::str("isogeny-perf.case/v1")),
                ("bits", Json::str(bits)),
                ("ell", Json::str(ell)),
                ("algorithm", Json::str(method.name())),
                ("input", input.clone()),
                ("input_sha256", Json::str(&input_sha)),
                ("output_sha256", Json::str(sha256_hex(expected.as_bytes()))),
                (
                    "count_model",
                    Json::str(
                        "generic Field calls; inversions opaque; setup and verification excluded",
                    ),
                ),
                ("operations", counts),
                ("diagnostic_ns", Json::Arr(raw)),
                ("qualified_timing", Json::Bool(false)),
                ("verified", Json::Bool(true)),
            ]);
            writeln!(out, "{}", row.dump()).map_err(|e| e.to_string())?;
            out.flush().map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}
fn main_run() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.first().is_some_and(|a| a == "compare") {
        return compare(&args[1..]);
    }
    if args.first().is_some_and(|a| a == "campaign") {
        return campaign::run(&args[1..]);
    }
    if args.is_empty() || args.iter().any(|a| a == "--help") {
        println!("isogeny-perf --out FRESH.jsonl --revision SOURCE_SHA [--bits 192,224,256,511] [--ells 7,31,101] [--samples 5] [--diagnostic-timing]\nisogeny-perf compare BASELINE.jsonl CANDIDATE.jsonl FRESH_COMPARISON.json\nAll eight BMSS methods; exact kernel/map equality and four point transports.\nCounts are an abstract arithmetic model. Timing is unreserved, diagnostic only.\nSynthetic supersingular fixtures at the requested field widths, not named NIST curves.");
        return Ok(());
    }
    let (mut path, mut revision, mut bits, mut ells, mut samples, mut timing) = (
        None,
        None,
        vec![192usize, 224],
        vec![7u64, 31],
        5usize,
        false,
    );
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if a == "--diagnostic-timing" {
            timing = true;
            i += 1;
            continue;
        }
        let value = args
            .get(i + 1)
            .ok_or_else(|| format!("missing value for {a}"))?;
        match a.as_str() {
            "--out" => path = Some(value.clone()),
            "--revision" => revision = Some(value.clone()),
            "--bits" => bits = csv(value)?,
            "--ells" => ells = csv(value)?,
            "--samples" => samples = value.parse().map_err(|_| "invalid samples")?,
            _ => return Err(format!("unknown argument {a}")),
        }
        i += 2;
    }
    let revision = revision.ok_or("--revision is required")?;
    if revision.len() != 40 || !revision.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err("revision must be a full 40-character git SHA".into());
    }
    if samples == 0 || samples > 100 {
        return Err("samples must be 1..100".into());
    }
    if bits.iter().any(|b| ![192, 224, 256, 511].contains(b)) {
        return Err("bits must be 192,224,256,511".into());
    }
    if ells
        .iter()
        .any(|l| *l < 5 || *l > 257 || !isogeny_algos::field::is_prime(*l))
    {
        return Err("ells must be primes in 5..257".into());
    }
    if bits
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>()
        .len()
        != bits.len()
        || ells
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != ells.len()
    {
        return Err("duplicate cases".into());
    }
    // ss_workload needs enough bits for 4*3*product(ells) and a random cofactor.
    let prime_bits = ells
        .iter()
        .map(|l| 64 - l.leading_zeros() as usize)
        .sum::<usize>()
        + 4;
    if bits.iter().any(|b| *b <= prime_bits + 16) {
        return Err("degree list is too large for field width".into());
    }
    let exe = std::env::current_exe().map_err(|e| e.to_string())?;
    let binary = std::fs::read(exe).map_err(|e| e.to_string())?;
    let mut out = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path.ok_or("--out is required")?)
        .map_err(|e| e.to_string())?;
    let header = Json::obj(vec![
        ("schema", Json::str("isogeny-perf.run/v1")),
        ("source_revision_supplied_by_launcher", Json::str(revision)),
        ("binary_sha256", Json::str(sha256_hex(&binary))),
        ("host", host()),
        ("args", Json::Arr(args.iter().map(Json::str).collect())),
    ]);
    writeln!(out, "{}", header.dump()).map_err(|e| e.to_string())?;
    for width in bits {
        match width {
            192 => run::<3>(&mut out, width, &ells, samples, timing)?,
            511 => run::<8>(&mut out, width, &ells, samples, timing)?,
            _ => run::<4>(&mut out, width, &ells, samples, timing)?,
        }
    }
    writeln!(
        out,
        "{}",
        Json::obj(vec![
            ("schema", Json::str("isogeny-perf.exit/v1")),
            ("status", Json::str("PASS"))
        ])
        .dump()
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}
fn text_field<'a>(j: &'a Json, key: &str) -> Result<&'a str, String> {
    j.get(key)
        .and_then(Json::as_str)
        .ok_or_else(|| format!("missing string field {key}"))
}
fn load_cases(path: &str) -> Result<(Json, Vec<Json>), String> {
    let raw = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut rows = raw
        .lines()
        .map(Json::parse)
        .collect::<Result<Vec<_>, _>>()?;
    if rows.len() < 3 {
        return Err("incomplete run".into());
    }
    let exit = rows.pop().unwrap();
    if text_field(&exit, "schema")? != "isogeny-perf.exit/v1"
        || text_field(&exit, "status")? != "PASS"
    {
        return Err("missing PASS exit receipt".into());
    }
    let header = rows.remove(0);
    if text_field(&header, "schema")? != "isogeny-perf.run/v1" {
        return Err("incorrect run schema".into());
    }
    let mut keys = std::collections::BTreeSet::new();
    for row in &rows {
        if text_field(row, "schema")? != "isogeny-perf.case/v1"
            || row.get("verified") != Some(&Json::Bool(true))
        {
            return Err("unverified case".into());
        }
        if !Method::ALL
            .iter()
            .any(|m| Some(m.name()) == row.get("algorithm").and_then(Json::as_str))
        {
            return Err("unknown algorithm".into());
        }
        let input = row.get("input").ok_or("missing input")?;
        if sha256_hex(input.dump().as_bytes()) != text_field(row, "input_sha256")? {
            return Err("input digest mismatch".into());
        }
        if !keys.insert((
            text_field(row, "bits")?,
            text_field(row, "ell")?,
            text_field(row, "algorithm")?,
        )) {
            return Err("duplicate case".into());
        }
        for op in ["add_sub", "mul", "inv", "neg"] {
            text_field(row.get("operations").ok_or("missing operations")?, op)?
                .parse::<u64>()
                .map_err(|_| "invalid operation count")?;
        }
    }
    // Each field-width/degree cell must contain all eight distinct BMSS methods.
    for row in &rows {
        for method in Method::ALL {
            if !keys.contains(&(
                text_field(row, "bits")?,
                text_field(row, "ell")?,
                method.name(),
            )) {
                return Err("missing algorithm in cell".into());
            }
        }
    }
    Ok((header, rows))
}
fn compare(args: &[String]) -> Result<(), String> {
    if args.len() != 3 {
        return Err("compare requires BASELINE CANDIDATE FRESH_OUTPUT".into());
    }
    let (ah, a) = load_cases(&args[0])?;
    let (bh, b) = load_cases(&args[1])?;
    if a.len() != b.len() {
        return Err("case count mismatch".into());
    }
    let mut pairs = vec![];
    for x in &a {
        let y = b
            .iter()
            .find(|y| {
                ["bits", "ell", "algorithm"]
                    .iter()
                    .all(|k| y.get(k) == x.get(k))
            })
            .ok_or("unmatched case")?;
        for key in ["input_sha256", "output_sha256", "count_model"] {
            if text_field(x, key)? != text_field(y, key)? {
                return Err(format!("{key} mismatch"));
            }
        }
        let mut ratios = vec![];
        for op in ["add_sub", "mul", "inv", "neg"] {
            let av = text_field(x.get("operations").unwrap(), op)?
                .parse::<u64>()
                .unwrap();
            let bv = text_field(y.get("operations").unwrap(), op)?
                .parse::<u64>()
                .unwrap();
            ratios.push((
                op,
                Json::obj(vec![
                    ("baseline", Json::str(av)),
                    ("candidate", Json::str(bv)),
                    (
                        "baseline_over_candidate",
                        if bv == 0 {
                            Json::Null
                        } else {
                            Json::Float(av as f64 / bv as f64)
                        },
                    ),
                ]),
            ));
        }
        pairs.push(Json::obj(vec![
            ("bits", x.get("bits").unwrap().clone()),
            ("ell", x.get("ell").unwrap().clone()),
            ("algorithm", x.get("algorithm").unwrap().clone()),
            ("input_sha256", x.get("input_sha256").unwrap().clone()),
            ("output_sha256", x.get("output_sha256").unwrap().clone()),
            ("operations", Json::obj(ratios)),
            ("qualified_wall_speedup", Json::Null),
        ]));
    }
    let summary = Json::obj(vec![
        ("schema", Json::str("isogeny-perf.comparison/v1")),
        ("status", Json::str("PASS")),
        (
            "baseline_run_sha256",
            Json::str(sha256_hex(
                &std::fs::read(&args[0]).map_err(|e| e.to_string())?,
            )),
        ),
        (
            "candidate_run_sha256",
            Json::str(sha256_hex(
                &std::fs::read(&args[1]).map_err(|e| e.to_string())?,
            )),
        ),
        (
            "baseline_host",
            ah.get("host").ok_or("missing host")?.clone(),
        ),
        (
            "candidate_host",
            bh.get("host").ok_or("missing host")?.clone(),
        ),
        ("exact_output_matches", Json::Num(pairs.len() as i64)),
        ("cases", Json::Arr(pairs)),
    ]);
    let mut f = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[2])
        .map_err(|e| e.to_string())?;
    writeln!(f, "{}", summary.dump()).map_err(|e| e.to_string())?;
    Ok(())
}
fn main() {
    if let Err(e) = main_run() {
        eprintln!("isogeny-perf: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn counts_are_separate_and_share_clones() {
        let f = Counted::new(isogeny_algos::field::Zp::new(101));
        assert_eq!(f.mul(3, 7), 21);
        assert_eq!(f.clone().div(21, 3), 7);
        assert_eq!(f.snapshot().get("mul").and_then(Json::as_str), Some("2"));
        assert_eq!(f.snapshot().get("inv").and_then(Json::as_str), Some("1"));
    }
}

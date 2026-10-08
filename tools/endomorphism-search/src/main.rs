use endomorphism_search::{
    curves::{self, Curve},
    orders::{self, Limits},
    Result, SCHEMA,
};
use serde_json::{json, Value};
use std::{collections::BTreeMap, env, fs, path::Path, process, time::Instant};

fn usage() -> &'static str {
    "endomorphism-search scan --discriminant=-619 [--degree-bound 1000] [--out FILE]\nendomorphism-search sweep [--max-abs-discriminant 1000] [--max-orders 5000]\nendomorphism-search probe --p 167 --a 25 --b 36 [--count-ops]\nendomorphism-search demo [--count-ops] [--out FILE]\nendomorphism-search replay --historical-dir DIR\nLimits: --max-work, --max-candidates, --seconds; probe also --max-kernels.\nNative counters are scalar stage diagnostics. Native timing and ECDLP speedup remain unset."
}
fn parse(args: &[String], allowed: &[&str]) -> Result<BTreeMap<String, String>> {
    let mut out = BTreeMap::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        let (key, value) = if let Some((k, v)) = a.split_once('=') {
            (k.to_string(), v.to_string())
        } else if a == "--count-ops" {
            (a.clone(), "true".into())
        } else {
            i += 1;
            (
                a.clone(),
                args.get(i).ok_or(format!("missing value for {a}"))?.clone(),
            )
        };
        if !allowed.contains(&key.as_str()) {
            return Err(format!("unknown option: {key}"));
        }
        if out.insert(key.clone(), value).is_some() {
            return Err(format!("duplicate option: {key}"));
        }
        i += 1;
    }
    Ok(out)
}
fn number<T: std::str::FromStr>(
    opts: &BTreeMap<String, String>,
    key: &str,
    default: Option<T>,
) -> Result<T> {
    match opts.get(key) {
        Some(v) => v.parse().map_err(|_| format!("invalid {key}: {v}")),
        None => default.ok_or(format!("required option {key}")),
    }
}
fn execute(args: &[String]) -> Result<(Value, Option<String>)> {
    let command = args.first().ok_or("a command is required")?;
    let allowed: Vec<_> = match command.as_str() {
        "scan" => vec![
            "--discriminant",
            "--degree-bound",
            "--max-work",
            "--max-candidates",
            "--seconds",
            "--out",
        ],
        "sweep" => vec![
            "--max-abs-discriminant",
            "--degree-bound",
            "--max-orders",
            "--max-work",
            "--max-candidates",
            "--seconds",
            "--out",
        ],
        "probe" => vec![
            "--p",
            "--a",
            "--b",
            "--degree-bound",
            "--max-kernels",
            "--seconds",
            "--count-ops",
            "--out",
        ],
        "demo" => vec!["--count-ops", "--out"],
        "replay" => vec!["--historical-dir", "--out"],
        _ => return Err(format!("unknown command: {command}")),
    };
    let opts = parse(&args[1..], &allowed)?;
    if opts.get("--count-ops").is_some_and(|v| v != "true") {
        return Err("--count-ops is a flag".into());
    }
    let bound = number(&opts, "--degree-bound", Some(1000i64))?;
    let limits = Limits {
        work: number(&opts, "--max-work", Some(1_000_000usize))?,
        candidates: number(&opts, "--max-candidates", Some(100_000usize))?,
        seconds: number(
            &opts,
            "--seconds",
            Some(if command == "sweep" { 30.0 } else { 10.0 }),
        )?,
    };
    limits.check()?;
    if !(1..=1_000_000_000_000).contains(&bound) {
        return Err("degree bound must be in 1..10^12".into());
    }
    let value = match command.as_str() {
        "scan" => orders::scan(number(&opts, "--discriminant", None)?, bound, limits)?,
        "sweep" => {
            let max = number(&opts, "--max-abs-discriminant", Some(1000i64))?;
            let max_orders = number(&opts, "--max-orders", Some(5000usize))?;
            if !(3..=100_000).contains(&max) || !(1..=50_000).contains(&max_orders) {
                return Err("sweep range must be 3..100000 and max-orders 1..50000".into());
            }
            let start = Instant::now();
            let mut stop = None;
            let mut rows = Vec::new();
            for size in 3..=max {
                if ![0, 1].contains(&(-size).rem_euclid(4)) {
                    continue;
                }
                let remaining = limits.seconds - start.elapsed().as_secs_f64();
                if rows.len() >= max_orders || remaining <= 0.0 {
                    stop = Some(if remaining <= 0.0 {
                        "time_limit"
                    } else {
                        "order_limit"
                    });
                    break;
                }
                rows.push(orders::summary(&orders::scan(
                    -size,
                    bound,
                    Limits {
                        seconds: remaining.min(2.0),
                        ..limits
                    },
                )?));
            }
            json!({"schema":SCHEMA,"discriminant_bound":max,"degree_bound":bound,"orders_scanned":rows.len(),"rows":rows,
                "status":if stop.is_none(){"complete_range_with_per_order_coverage"}else{"incomplete_range"},"stop_reason":stop,
                "limits":{"seconds":limits.seconds,"max_orders":max_orders},"wall_time":null})
        }
        "probe" => curves::probe(
            Curve::new(
                number(&opts, "--p", None)?,
                number(&opts, "--a", None)?,
                number(&opts, "--b", None)?,
            )?,
            bound,
            number(&opts, "--max-kernels", Some(32usize))?,
            limits.seconds,
            opts.contains_key("--count-ops"),
        )?,
        "demo" => {
            let mut results = Vec::new();
            for (p, a, b) in [(4093, 0, 2), (16381, 2, 0), (109, 43, 65), (167, 25, 36)] {
                results.push(curves::probe(
                    Curve::new(p, a, b)?,
                    1000,
                    32,
                    5.0,
                    opts.contains_key("--count-ops"),
                )?);
            }
            json!({"schema":SCHEMA,"curves":results,"scope":"synthetic small prime-field laboratory correctness; no native timing claim",
                "rho_speedup":null,"index_calculus_speedup":null})
        }
        "replay" => replay(Path::new(
            opts.get("--historical-dir")
                .ok_or("required --historical-dir")?,
        ))?,
        _ => unreachable!(),
    };
    Ok((value, opts.get("--out").cloned()))
}
fn read_json(path: &Path) -> Result<Value> {
    serde_json::from_slice(&fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?)
        .map_err(|e| e.to_string())
}
fn required_i(v: &Value, key: &str) -> Result<i64> {
    v.get(key)
        .and_then(Value::as_i64)
        .ok_or(format!("missing integer {key}"))
}
fn replay(dir: &Path) -> Result<Value> {
    // Pin the manifest itself: changing both a file and its claimed hash fails.
    // Correctness parity never attests historical or native timing.
    let manifest_bytes =
        fs::read(dir.join("results/artifact-manifest.json")).map_err(|e| e.to_string())?;
    if endomorphism_search::digest(&manifest_bytes)
        != "ca613e7ea9944ed4a93b9dfc11d5913dab3adf727fd1c2b0a2d02c6ad2cb4ee5"
    {
        return Err("historical manifest digest mismatch".into());
    }
    let manifest: Value = serde_json::from_slice(&manifest_bytes).map_err(|e| e.to_string())?;
    let entries = manifest["files"]
        .as_object()
        .ok_or("invalid artifact manifest")?;
    let mut hashes = 0;
    for (name, entry) in entries {
        if Path::new(name)
            .components()
            .any(|p| !matches!(p, std::path::Component::Normal(_)))
        {
            return Err("unsafe manifest path".into());
        }
        let bytes = fs::read(dir.join(name)).map_err(|e| e.to_string())?;
        if endomorphism_search::digest(&bytes) != entry.as_str().ok_or("manifest hash missing")? {
            return Err(format!("historical hash mismatch: {name}"));
        }
        hashes += 1;
    }
    let sweep = read_json(&dir.join("results/order-sweep.json"))?;
    let rows = sweep["rows"]
        .as_array()
        .ok_or("missing archived sweep rows")?;
    let mut orders_checked = 0;
    for row in rows {
        let native = orders::summary(&orders::scan(
            required_i(row, "discriminant")?,
            1000,
            Limits::default(),
        )?);
        for key in [
            "minimum_non_scalar_degree",
            "geometric_unit_count",
            "class_number",
            "abstract_candidate_count",
            "norm_search_status",
        ] {
            if native[key] != row[key] {
                return Err(format!(
                    "order parity mismatch: {} {key}",
                    row["discriminant"]
                ));
            }
        }
        orders_checked += 1;
    }
    if orders_checked != 500 {
        return Err("frozen sweep must contain all 500 orders".into());
    }
    let demo = read_json(&dir.join("results/demo.json"))?;
    let mut maps_checked = 0;
    let mut scalar_checks = 0;
    let curves = demo["curves"].as_array().ok_or("missing archived curves")?;
    if curves.len() != 4 {
        return Err("frozen demo must contain four curves".into());
    }
    for old in curves {
        let field = required_i(&old["curve"]["field"], "characteristic")?;
        let c = Curve::new(
            field,
            required_i(&old["curve"], "a")?,
            required_i(&old["curve"], "b")?,
        )?;
        let native = curves::probe(c, 1000, 32, 30.0, true)?;
        for key in ["point_count", "trace", "ordinary", "subgroup"] {
            if native[key] != old[key] {
                return Err(format!("curve parity mismatch: {key}"));
            }
        }
        if native["ring_binding"]["full_ring_discriminant"]
            != old["ring_binding"]["full_ring_discriminant"]
        {
            return Err("ring binding mismatch".into());
        }
        let native_maps = native["explicit_map_search"]["maps"]
            .as_array()
            .ok_or("missing native maps")?;
        let old_maps = old["explicit_map_search"]["maps"]
            .as_array()
            .ok_or("missing historical maps")?;
        if native_maps.len() != old_maps.len() {
            return Err("map count mismatch".into());
        }
        for map in old_maps {
            let new = native_maps
                .iter()
                .find(|m| m["map_id"] == map["map_id"])
                .ok_or("missing map")?;
            for key in [
                "family",
                "geometric_degree",
                "subgroup_eigenvalue",
                "trace_candidates",
            ] {
                if new[key] != map[key] {
                    return Err(format!("map parity mismatch: {} {key}", map["map_id"]));
                }
            }
            maps_checked += 1;
        }
        for map in native["scalar_operation_diagnostic"]["maps"]
            .as_array()
            .ok_or("missing scalar diagnostics")?
        {
            scalar_checks += map["exhaustive_scalar_checks"].as_u64().unwrap_or(0);
        }
    }
    Ok(
        json!({"schema":"endomorphism-search-replay/v1","historical_manifest_files_checked":hashes,
        "orders_checked":orders_checked,"curves_checked":curves.len(),"maps_checked":maps_checked,
        "exhaustive_accelerated_scalar_checks":scalar_checks,"status":"all_correctness_checks_passed",
        "historical_timings":"Python; preserved, not certified as native measurements","native_wall_timing":null,
        "speedup":null,"end_to_end_cost":null,
        "native_source_sha256":{
            "main.rs":endomorphism_search::digest(include_bytes!("main.rs")),
            "lib.rs":endomorphism_search::digest(include_bytes!("lib.rs")),
            "orders.rs":endomorphism_search::digest(include_bytes!("orders.rs")),
            "curves.rs":endomorphism_search::digest(include_bytes!("curves.rs")),
            "scalar.rs":endomorphism_search::digest(include_bytes!("scalar.rs"))}}),
    )
}
fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    if args == ["--help"] || args == ["-h"] {
        println!("{}", usage());
        return;
    }
    match execute(&args) {
        Ok((result, out)) => {
            let text = serde_json::to_string_pretty(&result).unwrap() + "\n";
            if let Some(out) = out {
                let path = Path::new(&out);
                let written = (|| -> std::io::Result<()> {
                    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
                        fs::create_dir_all(parent)?;
                    }
                    fs::write(path, text)
                })();
                if let Err(e) = written {
                    eprintln!("{out}: {e}");
                    process::exit(2);
                }
            } else {
                print!("{text}");
            }
        }
        Err(e) => {
            eprintln!("{e}\n{}", usage());
            process::exit(2);
        }
    }
}

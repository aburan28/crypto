//! Bounded, process-isolated construction search; see the frozen research protocol.
use isogeny_algos::{
    api::PrimeCurve, bigint::Big, field::is_prime, int::Int, json::Json, sha256::sha256_hex,
};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    thread,
    time::{Duration, Instant},
};

fn order(name: &str) -> Int {
    Int::from_big(&Big::from_dec(match name {
        "p192" => "6277101735386680763835789423176059013767194773182842284081",
        "p224" => "26959946667150639794667015087019625940457807714424391721682722368061",
        _ => panic!("unsupported preset"),
    }))
}
fn pow(mut a: u64, mut n: u64, p: u64) -> u64 {
    let mut r = 1;
    while n > 0 {
        if n & 1 != 0 {
            r = r * a % p;
        }
        a = a * a % p;
        n >>= 1;
    }
    r
}
fn mult_order(a: u64, p: u64) -> u64 {
    let mut r = 1;
    for n in 1..p {
        r = r * a % p;
        if r == 1 {
            return n;
        }
    }
    0
}
fn sieve(name: &str, upper: u64) -> Vec<(u64, i64, Json)> {
    let c = PrimeCurve::preset(name).unwrap();
    let t = &(&c.p + &Int::one()) - &order(name);
    (67..=upper)
        .filter(|&l| is_prime(l))
        .map(|l| {
            let tm = t.mod_u64(l);
            let pm = c.p.mod_u64(l);
            let d = (tm * tm + l - (4 * pm) % l) % l;
            let symbol = if d == 0 {
                0
            } else if pow(d, (l - 1) / 2, l) == 1 {
                1
            } else {
                -1
            };
            let roots: Vec<u64> = (0..l)
                .filter(|&x| (x * x + l - tm * x % l + pm) % l == 0)
                .collect();
            let row = Json::obj(vec![
                ("ell", Json::Num(l as i64)),
                ("discriminant_mod_ell", Json::Num(d as i64)),
                ("legendre", Json::Num(symbol)),
                (
                    "status",
                    Json::str(match symbol {
                        1 => "split",
                        -1 => "inert",
                        _ => "repeated_eigenvalue",
                    }),
                ),
                (
                    "eigenvalues",
                    Json::Arr(roots.iter().map(|&r| Json::Num(r as i64)).collect()),
                ),
                (
                    "eigenvalue_orders",
                    Json::Arr(
                        roots
                            .iter()
                            .map(|&r| Json::Num(mult_order(r, l) as i64))
                            .collect(),
                    ),
                ),
            ]);
            (l, symbol, row)
        })
        .collect()
}
fn write_json(path: &Path, value: &Json) {
    fs::write(path, format!("{}\n", value.dump())).unwrap();
}
fn run(cli: &Path, dir: &Path, args: &[String], stem: &str, timeout: u64) -> Json {
    let stdout_path = dir.join(format!("{stem}.json"));
    let stderr_path = dir.join(format!("{stem}.stderr.txt"));
    let mut child = Command::new(cli)
        .args(args)
        .stdout(Stdio::from(fs::File::create(&stdout_path).unwrap()))
        .stderr(Stdio::from(fs::File::create(&stderr_path).unwrap()))
        .spawn()
        .unwrap();
    let start = Instant::now();
    let mut timed_out = false;
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status;
        }
        if start.elapsed() >= Duration::from_secs(timeout) {
            timed_out = true;
            child.kill().unwrap();
            break child.wait().unwrap();
        }
        thread::sleep(Duration::from_millis(25));
    };
    let raw = fs::read_to_string(&stdout_path).unwrap();
    let parsed = Json::parse(&raw).ok();
    let n = parsed
        .as_ref()
        .and_then(|r| r.get("result"))
        .and_then(|r| r.get("isogenies"))
        .and_then(|r| {
            if let Json::Arr(v) = r {
                Some(v.len())
            } else {
                None
            }
        });
    let pass = !timed_out
        && status.success()
        && parsed
            .as_ref()
            .and_then(|r| r.get("status"))
            .and_then(Json::as_str)
            == Some("PASS");
    let construction = args.first().is_some_and(|a| a == "isogenies");
    let accepted = pass && (!construction || n == Some(2));
    let result = Json::obj(vec![
        (
            "command",
            Json::Arr(
                std::iter::once(cli.to_string_lossy().into_owned())
                    .chain(args.iter().cloned())
                    .map(Json::str)
                    .collect(),
            ),
        ),
        (
            "status",
            Json::str(if timed_out {
                "TIMEOUT"
            } else if accepted {
                "PASS"
            } else {
                "FAIL"
            }),
        ),
        (
            "exit_code",
            status.code().map_or(Json::Null, |c| Json::Num(c as i64)),
        ),
        (
            "elapsed_ms_operational_only",
            Json::Num(start.elapsed().as_millis() as i64),
        ),
        ("map_count", n.map_or(Json::Null, |n| Json::Num(n as i64))),
        (
            "stdout",
            Json::str(stdout_path.file_name().unwrap().to_string_lossy()),
        ),
        ("stdout_sha256", Json::str(sha256_hex(raw.as_bytes()))),
        (
            "stderr_sha256",
            Json::str(sha256_hex(&fs::read(&stderr_path).unwrap())),
        ),
    ]);
    write_json(&dir.join(format!("{stem}.receipt.json")), &result);
    println!(
        "{stem}: {} ({} ms)",
        result.get("status").unwrap().as_str().unwrap(),
        start.elapsed().as_millis()
    );
    std::io::stdout().flush().unwrap();
    result
}
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    assert!(
        args.len() == 3,
        "usage: search CLI OUTPUT_DIR SOURCE_COMMIT"
    );
    let cli = fs::canonicalize(&args[0]).unwrap();
    let out = PathBuf::from(&args[1]);
    fs::create_dir(&out).expect("a new evidence directory is required");
    let mut records = vec![];
    for name in ["p192", "p224"] {
        let rows = sieve(name, 4093);
        let mut degrees: Vec<u64> = rows
            .iter()
            .filter(|(l, s, _)| *l <= 257 && *s == 1)
            .map(|(l, _, _)| *l)
            .collect();
        for threshold in [509, 1009, 2003, 4000] {
            if let Some((l, _, _)) = rows.iter().find(|(l, s, _)| *l >= threshold && *s == 1) {
                degrees.push(*l);
            }
        }
        write_json(
            &out.join(format!("{name}-screen.json")),
            &Json::obj(vec![
                ("curve", Json::str(name)),
                ("source_order", Json::str(order(name))),
                (
                    "degrees_attempted",
                    Json::Arr(degrees.iter().map(|&l| Json::Num(l as i64)).collect()),
                ),
                (
                    "screen",
                    Json::Arr(rows.into_iter().map(|(_, _, r)| r).collect()),
                ),
            ]),
        );
        // Freeze the complete plan before starting this curve's construction attempts.
        let dir = out.join(name);
        fs::create_dir(&dir).unwrap();
        let receipt = run(
            &cli,
            &dir,
            &[
                "count".into(),
                "--curve".into(),
                name.into(),
                "--order".into(),
                order(name).to_string(),
                "--seed".into(),
                "1".into(),
            ],
            "order",
            180,
        );
        assert_eq!(
            receipt.get("status").unwrap().as_str(),
            Some("PASS"),
            "source order certification failed"
        );
        for ell in degrees {
            let receipt = run(
                &cli,
                &dir,
                &[
                    "isogenies".into(),
                    "--curve".into(),
                    name.into(),
                    "--ell".into(),
                    ell.to_string(),
                    "--order".into(),
                    order(name).to_string(),
                    "--seed".into(),
                    "1".into(),
                ],
                &format!("ell-{ell}"),
                180,
            );
            records.push(Json::obj(vec![
                ("curve", Json::str(name)),
                ("ell", Json::Num(ell as i64)),
                ("receipt", receipt),
            ]));
        }
    }
    write_json(
        &out.join("search.json"),
        &Json::obj(vec![
            ("schema", Json::str("large-degree-isogeny-search/v1")),
            ("source_commit", Json::str(&args[2])),
            (
                "protocol_sha256",
                Json::str(sha256_hex(include_bytes!(
                    "../../research/p192_p224_large_degree_isogenies_20261008/PROTOCOL.md"
                ))),
            ),
            ("timeout_seconds", Json::Num(180)),
            ("attempts", Json::Arr(records)),
        ]),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn split_and_inert_predictions_match_existing_small_degree_records() {
        // Direct roots and Euler's criterion must agree for every screened prime.
        for name in ["p192", "p224"] {
            for (_, s, row) in sieve(name, 4093) {
                let Json::Arr(roots) = row.get("eigenvalues").unwrap() else {
                    panic!()
                };
                assert_eq!(
                    roots.len(),
                    match s {
                        1 => 2,
                        -1 => 0,
                        _ => 1,
                    }
                );
            }
        }
    }
    #[test]
    fn published_orders_give_registry_traces() {
        for (name, trace) in [
            ("p192", "31607402316713927207482677199"),
            ("p224", "4733100108545601916421827343930821"),
        ] {
            let c = PrimeCurve::preset(name).unwrap();
            assert_eq!((&(&c.p + &Int::one()) - &order(name)).to_string(), trace);
        }
    }
}

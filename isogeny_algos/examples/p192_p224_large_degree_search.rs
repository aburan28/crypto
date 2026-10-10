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
                .filter(|&x| (x * x + l - tm * x % l + pm).is_multiple_of(l))
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
#[cfg(target_os = "macos")]
fn resident_bytes(pid: u32) -> Option<u64> {
    #[link(name = "proc")]
    extern "C" {
        fn proc_pidinfo(
            pid: i32,
            flavor: i32,
            arg: u64,
            buffer: *mut std::ffi::c_void,
            size: i32,
        ) -> i32;
    }
    let mut info = [0u64; 12];
    let size = std::mem::size_of_val(&info) as i32;
    // SDK proc_taskinfo: 96 bytes; the second u64 is resident size.
    let got = unsafe { proc_pidinfo(pid as i32, 4, 0, info.as_mut_ptr().cast(), size) };
    (got == size).then_some(info[1])
}
#[cfg(target_os = "linux")]
fn resident_bytes(pid: u32) -> Option<u64> {
    fs::read_to_string(format!("/proc/{pid}/status"))
        .ok()?
        .lines()
        .find_map(|line| {
            line.strip_prefix("VmRSS:")
                .and_then(|s| s.split_whitespace().next())
                .and_then(|s| s.parse::<u64>().ok())
                .map(|kb| kb * 1024)
        })
}
#[cfg(not(any(target_os = "macos", target_os = "linux")))]
fn resident_bytes(_: u32) -> Option<u64> {
    None
}
fn run(cli: &Path, dir: &Path, args: &[String], stem: &str, timeout: u64) -> Json {
    let stdout_path = dir.join(format!("{stem}.json"));
    let stderr_path = dir.join(format!("{stem}.stderr.txt"));
    let receipt_path = dir.join(format!("{stem}.receipt.json"));
    if receipt_path.exists() {
        let receipt = Json::parse(&fs::read_to_string(&receipt_path).unwrap()).unwrap();
        let expected_command = Json::Arr(
            std::iter::once(cli.to_string_lossy().into_owned())
                .chain(args.iter().cloned())
                .map(Json::str)
                .collect(),
        );
        assert_eq!(receipt.get("command"), Some(&expected_command));
        assert_eq!(
            receipt.get("stdout_sha256").and_then(Json::as_str),
            Some(sha256_hex(&fs::read(&stdout_path).unwrap()).as_str())
        );
        assert_eq!(
            receipt.get("stderr_sha256").and_then(Json::as_str),
            Some(sha256_hex(&fs::read(&stderr_path).unwrap()).as_str())
        );
        println!(
            "{stem}: reused sealed {} receipt",
            receipt.get("status").unwrap().as_str().unwrap()
        );
        return receipt;
    }
    if stdout_path.exists() || stderr_path.exists() {
        let interrupted = dir.join("interrupted");
        fs::create_dir_all(&interrupted).unwrap();
        let index = (1..)
            .find(|i| !interrupted.join(format!("{stem}-{i}")).exists())
            .unwrap();
        let archived = interrupted.join(format!("{stem}-{index}"));
        fs::create_dir(&archived).unwrap();
        let mut files = vec![];
        for path in [&stdout_path, &stderr_path] {
            if path.exists() {
                let bytes = fs::read(path).unwrap();
                files.push(Json::obj(vec![
                    (
                        "path",
                        Json::str(path.file_name().unwrap().to_string_lossy()),
                    ),
                    ("bytes", Json::Num(bytes.len() as i64)),
                    ("sha256", Json::str(sha256_hex(&bytes))),
                ]));
                fs::rename(path, archived.join(path.file_name().unwrap())).unwrap();
            }
        }
        write_json(&archived.join("interruption.json"), &Json::obj(vec![
            ("status", Json::str("INTERRUPTED_WITHOUT_RECEIPT")),
            ("elapsed_ms_operational_only", Json::Null),
            ("reason", Json::str("No surviving child process or completed receipt after session interruption; raw bytes preserved before rerun.")),
            ("files", Json::Arr(files)),
        ]));
        println!("{stem}: preserved interrupted output before retry");
    }
    let mut child = Command::new(cli)
        .args(args)
        .stdout(Stdio::from(fs::File::create(&stdout_path).unwrap()))
        .stderr(Stdio::from(fs::File::create(&stderr_path).unwrap()))
        .spawn()
        .unwrap();
    let start = Instant::now();
    let mut timed_out = false;
    let mut memory_limited = false;
    let mut monitor_unavailable = false;
    let mut peak_rss = 0;
    let mut unavailable_since = None;
    let mut unavailable_samples = 0;
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status;
        }
        if start.elapsed() >= Duration::from_secs(timeout) {
            timed_out = true;
            child.kill().unwrap();
            break child.wait().unwrap();
        }
        match resident_bytes(child.id()) {
            Some(rss) => {
                unavailable_since = None;
                peak_rss = peak_rss.max(rss);
                if rss > 8u64 << 30 {
                    memory_limited = true;
                    child.kill().unwrap();
                    break child.wait().unwrap();
                }
            }
            None => {
                if let Some(status) = child.try_wait().unwrap() {
                    break status;
                }
                unavailable_samples += 1;
                let first = unavailable_since.get_or_insert_with(Instant::now);
                // proc_taskinfo can disappear during exit before waitpid reports the
                // final status. Allow at most 100 ms for that transition; a live child
                // without restored monitoring still fails closed.
                if first.elapsed() >= Duration::from_millis(100) {
                    monitor_unavailable = true;
                    child.kill().unwrap();
                    break child.wait().unwrap();
                }
            }
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
        && !memory_limited
        && !monitor_unavailable
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
            Json::str(if memory_limited {
                "MEMORY_LIMIT"
            } else if monitor_unavailable {
                "RESOURCE_MONITOR_UNAVAILABLE"
            } else if timed_out {
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
        ("resident_memory_limit_bytes", Json::str(8u64 << 30)),
        ("sampled_peak_resident_bytes", Json::str(peak_rss)),
        ("memory_sampling_interval_ms", Json::Num(25)),
        ("unavailable_memory_samples", Json::Num(unavailable_samples)),
        ("memory_monitor_grace_ms", Json::Num(100)),
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
    write_json(&receipt_path, &result);
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
        args.len() == 3 || (args.len() == 4 && args[3] == "--resume"),
        "usage: search CLI OUTPUT_DIR SOURCE_COMMIT [--resume]"
    );
    let resume = args.len() == 4;
    let cli = fs::canonicalize(&args[0]).unwrap();
    let out = PathBuf::from(&args[1]);
    if resume {
        assert!(
            out.is_dir(),
            "resume requires an existing evidence directory"
        );
        assert!(
            !out.join("search.json").exists(),
            "completed search cannot be resumed"
        );
    } else {
        fs::create_dir(&out).expect("a new evidence directory is required");
    }
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
        let plan = Json::obj(vec![
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
        ]);
        let plan_path = out.join(format!("{name}-screen.json"));
        if plan_path.exists() {
            assert_eq!(
                fs::read_to_string(&plan_path).unwrap(),
                format!("{}\n", plan.dump()),
                "frozen search plan differs"
            );
        } else {
            write_json(&plan_path, &plan);
        }
        // Freeze the complete plan before starting this curve's construction attempts.
        let dir = out.join(name);
        fs::create_dir_all(&dir).unwrap();
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

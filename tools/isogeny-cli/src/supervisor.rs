//! Bounded, process-isolated construction search; see the frozen research protocol.
use isogeny_algos::{
    api::PrimeCurve, bigint::Big, field::is_prime, int::Int, json::Json, sha256::sha256_hex,
};
use std::{
    fs,
    io::Write,
    path::Path,
    process::{Command, Stdio},
    thread,
    time::{Duration, Instant},
};

pub fn order(name: &str) -> Int {
    Int::from_big(&Big::from_dec(match name {
        "p192" => "6277101735386680763835789423176059013767194773182842284081",
        "p224" => "26959946667150639794667015087019625940457807714424391721682722368061",
        _ => panic!("unsupported preset"),
    }))
}
pub fn pow(mut a: u64, mut n: u64, p: u64) -> u64 {
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
#[cfg(target_os = "windows")]
fn resident_bytes(pid: u32) -> Option<u64> {
    use std::ffi::c_void;
    #[repr(C)]
    struct Counters {
        cb: u32,
        faults: u32,
        peak: usize,
        working: usize,
        rest: [usize; 6],
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn OpenProcess(access: u32, inherit: i32, pid: u32) -> *mut c_void;
        fn CloseHandle(handle: *mut c_void) -> i32;
        fn K32GetProcessMemoryInfo(handle: *mut c_void, counters: *mut Counters, size: u32) -> i32;
    }
    // PROCESS_QUERY_LIMITED_INFORMATION. SDK PROCESS_MEMORY_COUNTERS layout.
    let mut info = Counters {
        cb: std::mem::size_of::<Counters>() as u32,
        faults: 0,
        peak: 0,
        working: 0,
        rest: [0; 6],
    };
    unsafe {
        let handle = OpenProcess(0x1000, 0, pid);
        if handle.is_null() {
            return None;
        }
        let ok = K32GetProcessMemoryInfo(handle, &mut info, std::mem::size_of::<Counters>() as u32);
        CloseHandle(handle);
        (ok != 0).then_some(info.working as u64)
    }
}
#[cfg(not(any(target_os = "macos", target_os = "linux", target_os = "windows")))]
fn resident_bytes(_: u32) -> Option<u64> {
    None
}
pub fn run(cli: &Path, dir: &Path, args: &[String], stem: &str, timeout: u64) -> Json {
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
    fs::write(dir.join(format!("{stem}.pid")), child.id().to_string()).unwrap();
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
    let bytes = fs::read(&stdout_path).unwrap();
    let raw = String::from_utf8_lossy(&bytes);
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
    let construction = args
        .first()
        .is_some_and(|a| a == "isogenies" || a == "kernel-first");
    let wanted = if args.first().is_some_and(|a| a == "kernel-first") {
        1
    } else {
        2
    };
    let accepted = pass && (!construction || n == Some(wanted));
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
        ("stdout_sha256", Json::str(sha256_hex(&bytes))),
        (
            "stderr_sha256",
            Json::str(sha256_hex(&fs::read(&stderr_path).unwrap())),
        ),
    ]);
    write_json(&receipt_path, &result);
    eprintln!(
        "{stem}: {} ({} ms)",
        result.get("status").unwrap().as_str().unwrap(),
        start.elapsed().as_millis()
    );
    std::io::stdout().flush().unwrap();
    result
}

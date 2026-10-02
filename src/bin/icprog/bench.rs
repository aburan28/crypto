//! The programme's runner (IC_TOOL_PROGRAM.md §5), native: the retired
//! `harness/bench.py`, step for step, so that a run tree it began can be
//! continued here without a seam.
//!
//! Every timed process is one `ic price --single-target` on one row, run
//! through the isolation tool (`isolated_bench run --wait --cpus 2`) with
//! `RAYON_NUM_THREADS=1`, after PSI `some avg10` has fallen below 4.0.  A
//! refused start is logged and retried after 15 s.  A contended or failed
//! process is kept and run again, at most twice; the first clean run is the
//! row's figure, and every attempt stays on disk.  An existing output is
//! never overwritten, so every command resumes where the last one stopped.
//!
//! Arms are named binaries.  `interleave` runs each row once per arm per
//! round, alternating the arms' order from round to round (A B, B A, ...).

use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use super::json::{self, J};
use super::runs;
use super::suite::Row;

pub const CPUS: &str = "2";
const REFUSAL_WAIT: Duration = Duration::from_secs(15);
const REFUSAL_LIMIT: usize = 60;
const PSI_READY: f64 = 4.0;
const PSI_WAIT_MAX: Duration = Duration::from_secs(120);

/// One arm: its name in the run tree and its binary.
pub struct Arm {
    pub name: String,
    pub binary: PathBuf,
}

/// The isolation tool every timed process runs under.
pub struct Bench {
    pub isolate: PathBuf,
}

/// The worst of the 10-second CPU and memory pressure (`some avg10`).
fn psi_some_avg10() -> f64 {
    let mut worst = 0.0f64;
    for kind in ["cpu", "memory"] {
        let Ok(text) = std::fs::read_to_string(format!("/proc/pressure/{kind}")) else {
            continue;
        };
        for line in text.lines().filter(|l| l.starts_with("some")) {
            if let Some(v) = line
                .split_once("avg10=")
                .and_then(|(_, rest)| rest.split_whitespace().next())
                .and_then(|v| v.parse::<f64>().ok())
            {
                worst = worst.max(v);
            }
        }
    }
    worst
}

/// Wait, at most two minutes, for the pressure to fall; the seconds waited.
fn wait_for_quiet() -> f64 {
    let start = Instant::now();
    while psi_some_avg10() >= PSI_READY && start.elapsed() < PSI_WAIT_MAX {
        std::thread::sleep(Duration::from_secs(1));
    }
    start.elapsed().as_secs_f64()
}

fn append(path: &Path, line: &str) -> Result<(), String> {
    let mut f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .map_err(|e| format!("{}: {e}", path.display()))?;
    f.write_all(line.as_bytes())
        .map_err(|e| format!("{}: {e}", path.display()))
}

/// UTC as `%Y-%m-%dT%H:%M:%SZ`.
fn utc_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs()) as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    // Civil date from days since 1970-01-01 (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

fn with_suffix_name(p: &Path, suffix: &str) -> PathBuf {
    let s = runs::stem(p);
    let name = format!("{}{suffix}", s.file_name().unwrap().to_string_lossy());
    s.with_file_name(name)
}

impl Bench {
    /// One timed process through the isolation tool, after PSI has fallen.
    pub fn launch(&self, cmd: &[String], out: &Path, log_dir: &Path) -> Result<(), String> {
        let rec = runs::record_path(out);
        let label = out
            .strip_prefix(log_dir)
            .map_err(|_| format!("{} is outside {}", out.display(), log_dir.display()))?
            .to_string_lossy()
            .into_owned();
        let err = with_suffix_name(out, ".stderr");
        let tmp = with_suffix_name(out, ".stderr.tmp");
        for attempt in 0..REFUSAL_LIMIT {
            let waited = wait_for_quiet();
            let stderr =
                std::fs::File::create(&tmp).map_err(|e| format!("{}: {e}", tmp.display()))?;
            Command::new(&self.isolate)
                .args(["run", "--wait", "--cpus", CPUS, "--out"])
                .arg(&rec)
                .args(["--label", &label, "--"])
                .args(cmd)
                .env("RAYON_NUM_THREADS", "1")
                .stdout(Stdio::null())
                .stderr(Stdio::from(stderr))
                .status()
                .map_err(|e| format!("{}: {e}", self.isolate.display()))?;
            if rec.exists() {
                std::fs::rename(&tmp, &err).map_err(|e| format!("{}: {e}", err.display()))?;
                append(
                    &log_dir.join("psi-waits.log"),
                    &format!("{label} {waited:.1}\n"),
                )?;
                return Ok(());
            }
            let why = std::fs::read_to_string(&tmp).unwrap_or_default();
            append(
                &log_dir.join("refusals.log"),
                &format!(
                    "{} {label} attempt {} (after {waited:.1} s of PSI wait): {}\n",
                    utc_now(),
                    attempt + 1,
                    why.trim()
                ),
            )?;
            let _ = std::fs::remove_file(&tmp);
            std::thread::sleep(REFUSAL_WAIT);
        }
        Err(format!(
            "{label}: the machine stayed busy for {REFUSAL_LIMIT} tries; stopped (resumable)"
        ))
    }

    /// The first clean run of up to three attempts; every attempt is kept.
    pub fn price(
        &self,
        binary: &Path,
        row: &Row,
        out: &Path,
        log_dir: &Path,
        extra: &[String],
    ) -> Result<J, String> {
        let mut rep = J::Null;
        for k in 0..=runs::RETRIES {
            let attempt = runs::attempt(out, k);
            if !attempt.exists() {
                if let Some(parent) = attempt.parent() {
                    std::fs::create_dir_all(parent)
                        .map_err(|e| format!("{}: {e}", parent.display()))?;
                }
                self.launch(&price_cmd(binary, row, &attempt, extra)?, &attempt, log_dir)?;
            }
            rep = runs::load(&attempt);
            if runs::clean(&attempt)? && rep.get("status").and_then(J::as_str) == Some("complete") {
                return Ok(rep);
            }
        }
        Ok(rep)
    }

    /// Every row, every round, every arm; the arms' order alternates by round.
    pub fn interleave(
        &self,
        arms: &[Arm],
        rows: &[Row],
        rounds: u32,
        out_dir: &Path,
        extra: &[String],
    ) -> Result<(), String> {
        for k in 1..=rounds {
            let order: Vec<&Arm> = if k % 2 == 1 {
                arms.iter().collect()
            } else {
                arms.iter().rev().collect()
            };
            for row in rows {
                for arm in &order {
                    let out = out_dir
                        .join(&arm.name)
                        .join(&row.id)
                        .join(format!("r{k}.price.json"));
                    let rep = self.price(&arm.binary, row, &out, out_dir, extra)?;
                    let status = rep
                        .get("status")
                        .and_then(J::as_str)
                        .unwrap_or("None")
                        .to_string();
                    let note = if status == "complete" {
                        let med = rep.get("median");
                        let f = |k: &str| med.and_then(|m| m.get(k)).and_then(J::as_f64);
                        format!(
                            "setup {:.3} online {:.4}",
                            f("s_setup").unwrap_or(f64::NAN),
                            f("s_ic_online").unwrap_or(f64::NAN)
                        )
                    } else {
                        String::new()
                    };
                    println!("r{k} {} {}: {status} {note}", row.id, arm.name);
                }
            }
        }
        Ok(())
    }
}

/// `ic price` on one row, as the suite's `command` gives it.
pub fn price_cmd(
    binary: &Path,
    row: &Row,
    out: &Path,
    extra: &[String],
) -> Result<Vec<String>, String> {
    let params = row
        .params
        .as_ref()
        .ok_or_else(|| format!("row {} has no parameter file", row.id))?;
    let rho_seed = row
        .rho_seed
        .ok_or_else(|| format!("row {} has no rho seed", row.id))?;
    let mut cmd = vec![
        binary.to_string_lossy().into_owned(),
        "price".into(),
        "--params".into(),
        params.to_string_lossy().into_owned(),
        "--json".into(),
        "--out".into(),
        out.to_string_lossy().into_owned(),
        "--single-target".into(),
        "--rho-seed".into(),
        rho_seed.to_string(),
    ];
    cmd.extend(extra.iter().cloned());
    Ok(cmd)
}

// ── the host manifest (AGENTS.md §10) ───────────────────────────────

fn sh(cmd: &str, args: &[&str]) -> String {
    Command::new(cmd)
        .args(args)
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

pub fn sha256_file(path: &Path) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(hex::encode(crypto_lib::hash::sha256::sha256(&bytes)))
}

fn glibc() -> Option<String> {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    {
        // SAFETY: gnu_get_libc_version returns a static NUL-terminated string.
        let v = unsafe { std::ffi::CStr::from_ptr(libc::gnu_get_libc_version()) };
        Some(v.to_string_lossy().into_owned())
    }
    #[cfg(not(all(target_os = "linux", target_env = "gnu")))]
    {
        None
    }
}

/// The host and binaries a run tree's processes ran on, written once.
pub fn host_manifest(
    out: &Path,
    root: &Path,
    binaries: J,
    isolate: &Path,
    note: &str,
) -> Result<J, String> {
    if out.exists() {
        return json::read(out);
    }
    let cpuinfo = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    let field = |key: &str| {
        cpuinfo
            .lines()
            .filter_map(|l| l.split_once(':'))
            .find(|(k, _)| k.trim() == key)
            .map(|(_, v)| v.trim().to_string())
    };
    let wanted: Vec<J> = field("flags")
        .unwrap_or_default()
        .split_whitespace()
        .filter(|f| {
            ["popcnt", "avx2", "pclmulqdq", "bmi2", "gfni", "vpclmulqdq"].contains(f)
                || f.starts_with("avx512")
        })
        .map(|f| J::Str(f.into()))
        .collect();
    let mem = std::fs::read_to_string("/proc/meminfo")
        .unwrap_or_default()
        .lines()
        .find(|l| l.starts_with("MemTotal"))
        .unwrap_or("")
        .to_string();
    let thp = |name: &str| {
        std::fs::read_to_string(format!("/sys/kernel/mm/transparent_hugepage/{name}"))
            .map_or(J::Null, |s| J::Str(s.trim().to_string()))
    };
    let cores = std::thread::available_parallelism().map_or(0, |n| n.get() as i128);
    let release = std::fs::read_to_string("/proc/sys/kernel/osrelease")
        .map(|s| s.trim().to_string())
        .unwrap_or_default();
    let arch = std::env::consts::ARCH;
    let libc_ver = glibc();
    let os = match &libc_ver {
        Some(v) => format!("Linux-{release}-{arch}-with-glibc{v}"),
        None => format!("Linux-{release}-{arch}"),
    };
    let lock = root.join("Cargo.lock");
    let siblings = (0..cores)
        .map(|c| {
            let p = format!("/sys/devices/system/cpu/cpu{c}/topology/thread_siblings_list");
            (
                c.to_string(),
                J::Str(
                    std::fs::read_to_string(p)
                        .unwrap_or_default()
                        .trim()
                        .to_string(),
                ),
            )
        })
        .collect();
    let root_s = root.to_string_lossy();
    let doc = J::Obj(vec![
        (
            "repository_commit".into(),
            J::Str(sh("git", &["-C", &root_s, "rev-parse", "HEAD"])),
        ),
        (
            "tree_status".into(),
            J::Str(sh("git", &["-C", &root_s, "status", "--porcelain"])),
        ),
        ("rustc".into(), J::Str(sh("rustc", &["--version"]))),
        (
            "cargo_lock_sha256".into(),
            if lock.exists() {
                J::Str(sha256_file(&lock)?)
            } else {
                J::Null
            },
        ),
        ("cpu_model".into(), field("model name").map_or(J::Null, J::Str)),
        ("cpu_flags_relevant".into(), J::Arr(wanted)),
        ("logical_cores".into(), J::Int(cores)),
        ("memory".into(), J::Str(mem)),
        (
            "transparent_hugepage".into(),
            J::Obj(vec![
                ("enabled".into(), thp("enabled")),
                ("defrag".into(), thp("defrag")),
            ]),
        ),
        ("os".into(), J::Str(os)),
        ("arch".into(), J::Str(arch.into())),
        (
            "glibc".into(),
            libc_ver.map_or(J::Null, |v| J::Str(format!("glibc {v}"))),
        ),
        ("binaries".into(), binaries),
        (
            "isolation_tool".into(),
            J::Obj(vec![
                ("path".into(), J::Str("src/bin/isolated_bench.rs".into())),
                ("binary_sha256".into(), J::Str(sha256_file(isolate)?)),
            ]),
        ),
        (
            "pinning".into(),
            J::Str(
                "isolated_bench run --wait --cpus 2 with RAYON_NUM_THREADS=1, after PSI some avg10 \
                 falls below 4.0; the tool's defaults (settle 2 s, other processes at most 0.10 \
                 CPUs, PSI at most 5); pins untimed under taskset -c 2"
                    .into(),
            ),
        ),
        ("smt_siblings".into(), J::Obj(siblings)),
        ("uptime".into(), J::Str(sh("uptime", &[]))),
        (
            "hardware_class".into(),
            J::Str("one x86-64 cloud container; no claim for Arm64, GPUs or other hosts".into()),
        ),
        ("runner".into(), J::Str("icprog (native)".into())),
        ("note".into(), J::Str(note.into())),
    ]);
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    std::fs::write(out, json::dumps(&doc, 1) + "\n")
        .map_err(|e| format!("{}: {e}", out.display()))?;
    Ok(doc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_utc_stamp_is_a_civil_date() {
        let s = utc_now();
        assert_eq!(s.len(), 20);
        assert!(s.ends_with('Z') && &s[4..5] == "-" && &s[10..11] == "T");
    }

    #[test]
    fn the_stderr_files_are_named_from_the_stem() {
        let out = Path::new("/x/cand/row/r3.price.json");
        assert_eq!(
            with_suffix_name(out, ".stderr"),
            Path::new("/x/cand/row/r3.stderr")
        );
        assert_eq!(
            with_suffix_name(out, ".stderr.tmp"),
            Path::new("/x/cand/row/r3.stderr.tmp")
        );
    }

    #[test]
    fn a_price_command_is_the_suites() {
        let row = Row {
            id: "s/M1-T01".into(),
            a: 1,
            n: 19,
            r: 262543.0,
            recipe_seed: Some(201),
            params: Some(PathBuf::from("/p/M1-T01.json")),
            rho_seed: Some(2293761),
        };
        let cmd = price_cmd(Path::new("/b/ic"), &row, Path::new("/o/r1.price.json"), &[]).unwrap();
        assert_eq!(
            cmd.join(" "),
            "/b/ic price --params /p/M1-T01.json --json --out /o/r1.price.json --single-target --rho-seed 2293761"
        );
    }
}

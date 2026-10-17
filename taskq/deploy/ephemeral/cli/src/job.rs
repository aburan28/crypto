//! Pure pieces of a job: ids, spans, quoting, the user data and the exit-status
//! mapping. Everything here must match what `../ephemeral` (bash) writes, since
//! `../runner.sh` reads the same files and either CLI can follow the other's jobs.

use anyhow::{Result, bail};
use serde::Serialize;

/// Phases after which a job never changes again.
pub const FINAL_PHASES: &[&str] = &[
    "done",
    "failed",
    "timeout",
    "interrupted",
    "error",
    "killed",
    "lost",
];

pub fn is_final(phase: &str) -> bool {
    FINAL_PHASES.contains(&phase)
}

/// Seconds in a span like 45m, 6h, 2d, 3600 or 3600s.
pub fn span_seconds(s: &str) -> Result<u64> {
    let (digits, unit) = match s.find(|c: char| !c.is_ascii_digit()) {
        Some(i) => s.split_at(i),
        None => (s, ""),
    };
    let mult = match unit {
        "" | "s" => 1,
        "m" => 60,
        "h" => 3600,
        "d" => 86400,
        _ => bail!("bad duration: {s} (use e.g. 90m, 6h, 2d)"),
    };
    if digits.is_empty() {
        bail!("bad duration: {s} (use e.g. 90m, 6h, 2d)");
    }
    let n: u64 = digits
        .parse()
        .map_err(|_| anyhow::anyhow!("bad duration: {s}"))?;
    Ok(n * mult)
}

/// The job-id label: `tr -cs 'A-Za-z0-9' '-' | tr A-Z a-z`, trimmed of dashes,
/// cut to 24 characters (the cut can leave a trailing dash, as in bash).
pub fn slug(name: &str) -> String {
    let mut out = String::new();
    let mut dash = false;
    for c in name.chars().chain(std::iter::once('\n')) {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
            dash = false;
        } else if !dash {
            out.push('-');
            dash = true;
        }
    }
    let trimmed = out.trim_matches('-');
    trimmed.chars().take(24).collect()
}

/// Default label: the basename of the command's first word.
pub fn default_name(first_arg: &str) -> String {
    let word = first_arg.split(' ').next().unwrap_or("");
    word.rsplit('/').next().unwrap_or(word).to_string()
}

/// Sortable, readable, unique: `20261017-041500-cargo-a1b2c3`.
pub fn job_id(utc_stamp: &str, name: &str, rand: [u8; 3]) -> String {
    let label = slug(name);
    let label = if label.is_empty() {
        "job".to_string()
    } else {
        label
    };
    format!(
        "{utc_stamp}-{label}-{:02x}{:02x}{:02x}",
        rand[0], rand[1], rand[2]
    )
}

pub fn valid_id(id: &str) -> bool {
    !id.is_empty()
        && id
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-'))
}

/// Quote one word for bash. Plain words stay as they are, like `printf %q`.
pub fn shell_quote(s: &str) -> String {
    if !s.is_empty()
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || "_@%+=:,./-".contains(c))
    {
        return s.to_string();
    }
    format!("'{}'", s.replace('\'', r"'\''"))
}

/// `cmd.sh`: one argument is a shell snippet, several are an argv.
pub fn cmd_script(args: &[String]) -> String {
    if let [snippet] = args {
        return format!("{snippet}\n");
    }
    let mut s = String::from("exec");
    for a in args {
        s.push(' ');
        s.push_str(&shell_quote(a));
    }
    s.push('\n');
    s
}

/// The `env` file the runner sources: `KEY=quoted-value` per line.
pub fn env_file(envs: &[String]) -> String {
    envs.iter()
        .filter_map(|kv| kv.split_once('='))
        .map(|(k, v)| format!("{k}={}\n", shell_quote(v)))
        .collect()
}

/// User data. Its first act arms the deadline, before anything else can fail.
pub fn user_data(id: &str, bucket: &str, region: &str, deadline: u64, keep: bool) -> String {
    format!(
        r#"#!/bin/bash
left=$(( {deadline} - $(date +%s) ))
[ "$left" -gt 0 ] || shutdown -h now
systemd-run --unit=ephemeral-deadline --on-active="${{left}}s" /sbin/shutdown -h now
export EPHEMERAL_JOB={id} EPHEMERAL_BUCKET={bucket} EPHEMERAL_DEADLINE={deadline} AWS_DEFAULT_REGION={region}
for i in $(seq 1 20); do
    aws s3 cp --only-show-errors s3://{bucket}/jobs/{id}/runner.sh /root/ephemeral-runner.sh && break
    sleep 6
done
bash /root/ephemeral-runner.sh
{keep} || shutdown -h now
"#
    )
}

/// Exit status of `run` and `logs -f` for a final phase.
pub fn exit_status(phase: &str, code: Option<i64>) -> i32 {
    match phase {
        "done" => 0,
        "failed" => code.map(|c| c as i32).unwrap_or(1),
        "timeout" => 124,
        "interrupted" => 75,
        _ => 70,
    }
}

/// https form of a GitHub ssh remote; anything else unchanged.
pub fn origin_https(url: &str) -> String {
    if let Some(rest) = url.strip_prefix("git@github.com:") {
        format!("https://github.com/{rest}")
    } else if let Some(rest) = url.strip_prefix("ssh://git@github.com/") {
        format!("https://github.com/{rest}")
    } else {
        url.to_string()
    }
}

/// `--cargo-dir` must stay inside the checkout.
pub fn check_cargo_dir(dir: &str) -> Result<String> {
    let dir = dir.trim_end_matches('/');
    let dir = if dir.is_empty() { "." } else { dir };
    if dir.starts_with('/') || format!("/{dir}/").contains("/../") {
        bail!("--cargo-dir is relative to the checkout");
    }
    Ok(dir.to_string())
}

pub fn human_size(n: i64) -> String {
    let units = ["Bytes", "KiB", "MiB", "GiB", "TiB"];
    let mut v = n as f64;
    let mut u = 0;
    while v >= 1024.0 && u + 1 < units.len() {
        v /= 1024.0;
        u += 1;
    }
    if u == 0 {
        format!("{n} Bytes")
    } else {
        format!("{v:.1} {}", units[u])
    }
}

/// `job.json`, field for field as the bash CLI writes it.
#[derive(Serialize, Debug, Clone)]
pub struct JobSpec {
    pub id: String,
    pub created: String,
    pub owner: String,
    pub instance_type: String,
    pub spot: bool,
    pub ttl: String,
    pub deadline: u64,
    pub disk_gb: i32,
    pub ami: String,
    pub toolchain: String,
    pub setup: String,
    pub repo_url: String,
    pub commit: String,
    pub local_head: String,
    pub patch: bool,
    pub keep: bool,
    pub command_display: String,
    pub paths: Vec<String>,
    pub artifacts: Vec<String>,
    pub cache: bool,
    pub cargo_dir: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spans() {
        assert_eq!(span_seconds("3600").unwrap(), 3600);
        assert_eq!(span_seconds("45s").unwrap(), 45);
        assert_eq!(span_seconds("90m").unwrap(), 5400);
        assert_eq!(span_seconds("6h").unwrap(), 21600);
        assert_eq!(span_seconds("2d").unwrap(), 172800);
        for bad in ["", "h", "6x", "1.5h", "-1h", "6hh"] {
            assert!(span_seconds(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn slugs_match_tr() {
        assert_eq!(slug("cargo"), "cargo");
        assert_eq!(slug("Cargo Run!!x"), "cargo-run-x");
        assert_eq!(slug("--a__b--"), "a-b");
        assert_eq!(slug(""), "");
        assert_eq!(
            slug("abcdefghijklmnopqrstuvw-yz"),
            "abcdefghijklmnopqrstuvw-"
        );
        assert_eq!(default_name("/usr/bin/time -v cargo"), "time");
        assert_eq!(
            job_id("20261017-041500", "", [0xab, 0x01, 0xff]),
            "20261017-041500-job-ab01ff"
        );
        assert_eq!(
            job_id("20261017-041500", "bake-x86_64", [1, 2, 3]),
            "20261017-041500-bake-x86-64-010203"
        );
        assert!(valid_id("20261017-041500-job-ab01ff"));
        assert!(!valid_id("../x"));
        assert!(!valid_id("a/b"));
        assert!(!valid_id(""));
    }

    #[test]
    fn quoting_round_trips_through_bash() {
        let args: Vec<String> = ["echo", "a b", "it's", "", "$HOME", "x;y"]
            .map(String::from)
            .to_vec();
        let script = cmd_script(&args);
        assert!(script.starts_with("exec echo 'a b' 'it'\\''s' '' '$HOME' 'x;y'"));
        let out = std::process::Command::new("bash")
            .arg("-c")
            .arg(
                script
                    .trim_start_matches("exec ")
                    .replace("echo", "printf '[%s]'"),
            )
            .output()
            .unwrap();
        assert_eq!(
            String::from_utf8_lossy(&out.stdout),
            "[a b][it's][][$HOME][x;y]"
        );
        assert_eq!(cmd_script(&["ls | wc -l".into()]), "ls | wc -l\n");
        assert_eq!(
            env_file(&["A=1".into(), "B=x y=z".into()]),
            "A=1\nB='x y=z'\n"
        );
    }

    #[test]
    fn user_data_arms_deadline_first() {
        let ud = user_data("j1", "bkt", "us-east-1", 1_800_000_000, false);
        let lines: Vec<&str> = ud.lines().collect();
        assert_eq!(lines[0], "#!/bin/bash");
        assert_eq!(lines[1], "left=$(( 1800000000 - $(date +%s) ))");
        assert!(lines[3].contains(r#"--on-active="${left}s""#));
        assert!(ud.contains("s3://bkt/jobs/j1/runner.sh"));
        assert!(ud.ends_with("false || shutdown -h now\n"));
        let out = std::process::Command::new("bash")
            .arg("-n")
            .arg("-c")
            .arg(&ud)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }

    #[test]
    fn exit_statuses() {
        assert_eq!(exit_status("done", Some(0)), 0);
        assert_eq!(exit_status("failed", Some(7)), 7);
        assert_eq!(exit_status("failed", None), 1);
        assert_eq!(exit_status("timeout", None), 124);
        assert_eq!(exit_status("interrupted", None), 75);
        assert_eq!(exit_status("lost", None), 70);
        assert_eq!(exit_status("error", None), 70);
    }

    #[test]
    fn misc() {
        assert_eq!(
            origin_https("git@github.com:aburan28/crypto.git"),
            "https://github.com/aburan28/crypto.git"
        );
        assert_eq!(
            origin_https("ssh://git@github.com/a/b"),
            "https://github.com/a/b"
        );
        assert_eq!(origin_https("https://x/y"), "https://x/y");
        assert_eq!(check_cargo_dir("tools/x/").unwrap(), "tools/x");
        assert_eq!(check_cargo_dir(".").unwrap(), ".");
        assert!(check_cargo_dir("/abs").is_err());
        assert!(check_cargo_dir("a/../../b").is_err());
        assert!(check_cargo_dir("..").is_err());
        assert_eq!(human_size(10), "10 Bytes");
        assert_eq!(human_size(1536), "1.5 KiB");
        assert_eq!(human_size(795 << 20), "795.0 MiB");
    }
}

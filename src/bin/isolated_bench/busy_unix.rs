//! Native serialization for build/test work on Unix. This is not a timing harness.
//! Can be bootstrapped with rustc before Cargo builds the full isolation tool.
use std::{fs::OpenOptions, os::fd::AsRawFd, path::Path, process::Command};

unsafe extern "C" {
    fn flock(fd: i32, operation: i32) -> i32;
}

pub fn run(lock: &Path, command: &[String]) -> Result<i32, String> {
    let (program, arguments) = command.split_first().ok_or("busy needs a command")?;
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(lock)
        .map_err(|e| format!("open {}: {e}", lock.display()))?;
    // LOCK_EX is 2 on macOS and Linux. The file remains open until the child exits.
    loop {
        if unsafe { flock(file.as_raw_fd(), 2) } == 0 {
            break;
        }
        let error = std::io::Error::last_os_error();
        if error.kind() != std::io::ErrorKind::Interrupted {
            return Err(format!("lock {}: {error}", lock.display()));
        }
    }
    let status = Command::new(program)
        .args(arguments)
        .status()
        .map_err(|e| format!("launch {program}: {e}"))?;
    Ok(status.code().unwrap_or(1))
}

pub fn cli() -> std::process::ExitCode {
    let mut args = std::env::args().skip(1).collect::<Vec<_>>();
    let mut lock = std::path::PathBuf::from("/tmp/crypto-bench.lock");
    if args.first().is_some_and(|a| a == "busy") {
        args.remove(0);
    } else {
        eprintln!("Only busy is available here. Timed CPU isolation requires Linux.");
        return std::process::ExitCode::FAILURE;
    }
    if args.first().is_some_and(|a| a == "--lock") && args.len() >= 2 {
        lock = args[1].clone().into();
        args.drain(..2);
    }
    if args.first().is_some_and(|a| a == "--") {
        args.remove(0);
    } else {
        eprintln!("usage: isolated_bench busy [--lock PATH] -- COMMAND [ARGS]");
        return std::process::ExitCode::FAILURE;
    }
    match run(&lock, &args) {
        Ok(code) => std::process::ExitCode::from(code.clamp(0, 255) as u8),
        Err(error) => {
            eprintln!("isolated_bench: {error}");
            std::process::ExitCode::FAILURE
        }
    }
}

// The standalone entrypoint is unused when included by isolated_bench.
#[allow(dead_code)]
fn main() -> std::process::ExitCode {
    cli()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

    struct Directory(std::path::PathBuf);
    impl Directory {
        fn new() -> Self {
            let name = format!(
                "ic-native-busy-{}-{}",
                std::process::id(),
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            );
            let path = std::env::temp_dir().join(name);
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for Directory {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn busy_serializes_two_children_on_the_same_lock() {
        let directory = Directory::new();
        let lock = directory.0.join("lock");
        let ready = directory.0.join("ready");
        let events = directory.0.join("events");
        let first_lock = lock.clone();
        let first_ready = ready.to_str().unwrap().to_string();
        let first_events = events.to_str().unwrap().to_string();
        let child = std::thread::spawn(move || {
            run(
                &first_lock,
                &[
                    "/bin/sh".into(),
                    "-c".into(),
                    "printf ready > \"$1\"; sleep 0.2; printf A >> \"$2\"".into(),
                    "busy-test".into(),
                    first_ready,
                    first_events,
                ],
            )
        });
        let deadline = Instant::now() + Duration::from_secs(10);
        while !ready.exists() && Instant::now() < deadline {
            std::thread::sleep(Duration::from_millis(5));
        }
        assert!(ready.exists(), "first child did not start");
        assert_eq!(
            run(
                &lock,
                &[
                    "/bin/sh".into(),
                    "-c".into(),
                    "printf B >> \"$1\"".into(),
                    "busy-test".into(),
                    events.to_str().unwrap().into()
                ]
            )
            .unwrap(),
            0
        );
        assert_eq!(child.join().unwrap().unwrap(), 0);
        assert_eq!(std::fs::read_to_string(&events).unwrap(), "AB");
    }

    #[test]
    fn status_and_launch_errors_release_the_lock() {
        let directory = Directory::new();
        let lock = directory.0.join("lock");
        assert_eq!(
            run(&lock, &["/bin/sh".into(), "-c".into(), "exit 3".into()]).unwrap(),
            3
        );
        assert!(run(
            &lock,
            &[directory.0.join("absent").to_str().unwrap().into()]
        )
        .is_err());
        assert_eq!(
            run(&lock, &["/bin/sh".into(), "-c".into(), "exit 0".into()]).unwrap(),
            0
        );
    }
}

//! Disclosed-input CryptoMiniSat prestarted-stdin control, not an IC run.
//!
//! Run only after the timed ordinary panels release the shared busy lock:
//!   cargo run --locked --release --example ic_sat_stdin_probe -- \
//!     CMS CNF ANF MANIFEST PREPARATION POINT_X POINT_Y NEW_OUTPUT_DIR \
//!     EXPECTED_CMS_SHA256 EXPECTED_MANIFEST_SHA256 CONFLICT_BUDGET TIMEOUT_MS \
//!     EXPECTED_STATUS
//!
//! The accepted CMS CLI receives the same frozen CNF through a regular file
//! and through its prestarted standard-input reader. This proves only transport compatibility on
//! one disclosed source instance; it does not admit a SAT target solver.

#[cfg(unix)]
mod unix_probe {
    use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
    use crypto_lib::cryptanalysis::prepared_sat_control::{
        native_status, parse_model, sha256, validate_manifest, validate_source_syntax, verify_anf,
        verify_cnf, NativeOutput, PreparedState, ORDER,
    };
    use serde_json::{json, Value};
    use std::{
        collections::HashMap,
        env, fs,
        fs::OpenOptions,
        io::{ErrorKind, Write},
        os::unix::{io::AsRawFd, process::CommandExt},
        path::{Path, PathBuf},
        process::{Child, ChildStdin, Command, ExitCode, Stdio},
        time::{Duration, Instant},
    };

    const MAX_CNF_BYTES: usize = 16 * 1024 * 1024;
    const MAX_ANF_BYTES: usize = 16 * 1024 * 1024;
    const MAX_OUTPUT_BYTES: usize = 8 * 1024 * 1024;

    fn read_bounded(path: &Path, limit: usize) -> Result<Vec<u8>, String> {
        let len = fs::metadata(path).map_err(|e| e.to_string())?.len();
        if len > limit as u64 {
            return Err(format!("{} exceeds {limit} bytes", path.display()));
        }
        fs::read(path).map_err(|e| e.to_string())
    }

    fn save_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(|e| e.to_string())?;
        file.write_all(bytes).map_err(|e| e.to_string())?;
        file.sync_all().map_err(|e| e.to_string())
    }

    struct ChildGuard(Child);
    impl Drop for ChildGuard {
        fn drop(&mut self) {
            // A solver can fork and exit while descendants still hold resources.
            // Reap its entire dedicated process group even after parent exit.
            let pid = self.0.id();
            if pid > 1 {
                unsafe { libc::kill(-(pid as i32), libc::SIGKILL) };
            }
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }

    fn args(input: Option<&Path>, conflict_budget: u64) -> Vec<String> {
        let mut args = vec![
            "--verb".into(),
            "1".into(),
            "--threads".into(),
            "1".into(),
            "--random".into(),
            "1".into(),
            "--maxsol".into(),
            "1".into(),
            "--maxconfl".into(),
            conflict_budget.to_string(),
        ];
        if let Some(input) = input {
            args.push(input.to_string_lossy().into_owned());
        }
        args
    }

    fn launch(
        cms: &Path,
        input: Option<&Path>,
        out: &Path,
        stem: &str,
        conflict_budget: u64,
    ) -> Result<(ChildGuard, Vec<String>), String> {
        let argv = args(input, conflict_budget);
        let stdout = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join(format!("{stem}.stdout")))
            .map_err(|e| e.to_string())?;
        let stderr = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(out.join(format!("{stem}.stderr")))
            .map_err(|e| e.to_string())?;
        let mut command = Command::new(cms);
        command
            .args(&argv)
            .current_dir(out)
            .env_clear()
            .env("LC_ALL", "C")
            .stdin(if input.is_some() {
                Stdio::null()
            } else {
                Stdio::piped()
            })
            .stdout(Stdio::from(stdout))
            .stderr(Stdio::from(stderr))
            .process_group(0);
        let child = command.spawn().map_err(|e| e.to_string())?;
        Ok((ChildGuard(child), argv))
    }

    fn wait(child: &mut ChildGuard, timeout: Duration) -> Result<(Option<i32>, bool), String> {
        let start = Instant::now();
        loop {
            if let Some(status) = child.0.try_wait().map_err(|e| e.to_string())? {
                return Ok((status.code(), false));
            }
            if start.elapsed() >= timeout {
                return Ok((None, true));
            }
            std::thread::sleep(Duration::from_millis(5));
        }
    }

    fn geometric_witness(state: &PreparedState, public: FastPoint) -> Option<[usize; 3]> {
        let mut pairs = HashMap::new();
        for (i, &a) in state.geometry.iter().enumerate() {
            for (j, &b) in state.geometry.iter().enumerate() {
                pairs.entry(state.curve.add(a, b)).or_insert((i, j));
            }
        }
        for (k, &c) in state.geometry.iter().enumerate() {
            let wanted = state.curve.add(public, state.curve.neg(c));
            if let Some(&(i, j)) = pairs.get(&wanted) {
                if state
                    .curve
                    .add(state.curve.add(state.geometry[i], state.geometry[j]), c)
                    == public
                {
                    return Some([i, j, k]);
                }
            }
        }
        None
    }

    fn model_check(
        status: &str,
        stdout: &str,
        anf: &str,
        cnf: &str,
        state: &PreparedState,
        public: FastPoint,
    ) -> Result<Value, String> {
        if status != "SAT_MODEL" {
            return Ok(Value::Null);
        }
        let width = cnf
            .lines()
            .find(|line| line.starts_with("p cnf "))
            .and_then(|line| line.split_whitespace().nth(2))
            .and_then(|value| value.parse::<usize>().ok())
            .ok_or("CNF variable count is missing")?;
        let model = parse_model(stdout, width)?;
        verify_anf(anf, model.get(..51).ok_or("short SAT model")?)?;
        verify_cnf(cnf, &model)?;
        let indices = state
            .lift(&model, public)
            .ok_or("source-valid SAT model does not lift to a full-point decomposition")?;
        Ok(json!({"source_model_valid":true,"model_sha256":sha256(
        &model.iter().map(|&v| u8::from(v)).collect::<Vec<_>>()),
        "full_point_witness_indices":indices}))
    }

    fn outcome(
        out: &Path,
        stem: &str,
        argv: &[String],
        pid: u32,
        result: (Option<i32>, bool),
        anf: &str,
        cnf: &str,
        state: &PreparedState,
        public: FastPoint,
    ) -> Result<Value, String> {
        let stdout_bytes = read_bounded(&out.join(format!("{stem}.stdout")), MAX_OUTPUT_BYTES)?;
        let stderr_bytes = read_bounded(&out.join(format!("{stem}.stderr")), MAX_OUTPUT_BYTES)?;
        let stdout = String::from_utf8(stdout_bytes.clone()).map_err(|e| e.to_string())?;
        let native = NativeOutput {
            exit_code: result.0,
            timed_out: result.1,
            stdout: stdout.clone(),
            receipt: Value::Null,
        };
        let status = native_status(&native);
        Ok(
            json!({"pid":pid,"argv":argv,"exit_code":result.0,"timed_out":result.1,
        "status":status,"model_check":model_check(status, &stdout, anf, cnf, state, public)?,
        "stdout_bytes":stdout_bytes.len(),"stdout_sha256":sha256(&stdout_bytes),
        "stderr_bytes":stderr_bytes.len(),"stderr_sha256":sha256(&stderr_bytes)}),
        )
    }

    fn reader_ready(out: &Path, child: &mut ChildGuard, deadline: Instant) -> Result<u64, String> {
        let start = Instant::now();
        let stdout_path = out.join("stdin.stdout");
        const READY: &str = "c Reading from standard input... Use '-h' or '--help' for help.";
        loop {
            let bytes = read_bounded(&stdout_path, MAX_OUTPUT_BYTES)?;
            if String::from_utf8_lossy(&bytes)
                .lines()
                .any(|line| line == READY)
            {
                return u64::try_from(start.elapsed().as_nanos())
                    .map_err(|_| "reader-ready duration overflow".into());
            }
            if let Some(status) = child.0.try_wait().map_err(|e| e.to_string())? {
                return Err(format!(
                    "CMS exited before stdin reader-ready line: {status}"
                ));
            }
            if Instant::now() >= deadline {
                return Err("CMS did not reach its stdin reader before deadline".into());
            }
            std::thread::sleep(Duration::from_millis(5));
        }
    }

    fn send_stdin(
        writer: &mut ChildStdin,
        bytes: &[u8],
        child: &mut ChildGuard,
        deadline: Instant,
    ) -> Result<(), String> {
        let fd = writer.as_raw_fd();
        let flags = unsafe { libc::fcntl(fd, libc::F_GETFL) };
        if flags < 0 || unsafe { libc::fcntl(fd, libc::F_SETFL, flags | libc::O_NONBLOCK) } < 0 {
            return Err("cannot enable nonblocking stdin delivery".into());
        }
        let mut sent = 0;
        while sent < bytes.len() {
            match writer.write(&bytes[sent..bytes.len().min(sent + 65_536)]) {
                Ok(0) => return Err("stdin writer made no progress".into()),
                Ok(n) => sent += n,
                Err(error)
                    if matches!(error.kind(), ErrorKind::WouldBlock | ErrorKind::Interrupted) => {}
                Err(error) => return Err(format!("stdin write: {error}")),
            }
            if child.0.try_wait().map_err(|e| e.to_string())?.is_some() && sent < bytes.len() {
                return Err("CMS exited before complete stdin CNF delivery".into());
            }
            if Instant::now() >= deadline {
                return Err("stdin delivery exceeded deadline".into());
            }
            if sent < bytes.len() {
                std::thread::sleep(Duration::from_millis(1));
            }
        }
        Ok(())
    }

    fn run(
        cms: &Path,
        cnf_path: &Path,
        anf_path: &Path,
        manifest_path: &Path,
        preparation_path: &Path,
        point: [u64; 2],
        out: &Path,
        expected: &str,
        expected_manifest: &str,
        conflict_budget: u64,
        timeout_ms: u64,
        expected_status: &str,
    ) -> Result<Value, String> {
        if ![
            cms,
            cnf_path,
            anf_path,
            manifest_path,
            preparation_path,
            out,
        ]
        .iter()
        .all(|path| path.is_absolute())
        {
            return Err("CMS, source and output paths must be absolute".into());
        }
        let source_dir = manifest_path.parent().ok_or("manifest lacks parent")?;
        if cnf_path.parent() != Some(source_dir)
            || anf_path.parent() != Some(source_dir)
            || cnf_path.file_name() != Some(std::ffi::OsStr::new("instance.xor.cnf"))
            || anf_path.file_name() != Some(std::ffi::OsStr::new("instance.anf"))
            || manifest_path.file_name() != Some(std::ffi::OsStr::new("manifest.json"))
        {
            return Err("source inputs must be one retained exporter instance".into());
        }
        for (name, value) in [("CMS", expected), ("manifest", expected_manifest)] {
            if value.len() != 64
                || !value
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            {
                return Err(format!(
                    "expected {name} SHA-256 must be 64 lowercase hex digits"
                ));
            }
        }
        if !(100..=120_000).contains(&timeout_ms) {
            return Err("timeout must be 100..120000 ms".into());
        }
        if !(1..=1_000_000).contains(&conflict_budget) {
            return Err("conflict budget must be 1..1000000".into());
        }
        if !matches!(
            expected_status,
            "SAT_MODEL" | "SOURCE_UNSAT" | "CONFLICT_BUDGET_INCONCLUSIVE" | "UNKNOWN_INCONCLUSIVE"
        ) {
            return Err("expected status must be a recognized CMS disposition".into());
        }
        let cms_bytes = read_bounded(cms, 8 * 1024 * 1024)?;
        let cms_sha = sha256(&cms_bytes);
        if cms_sha != expected {
            return Err("CMS executable SHA-256 differs from expected pin".into());
        }
        let cnf_bytes = read_bounded(cnf_path, MAX_CNF_BYTES)?;
        let cnf = String::from_utf8(cnf_bytes.clone()).map_err(|e| e.to_string())?;
        let anf_bytes = read_bounded(anf_path, MAX_ANF_BYTES)?;
        let anf = String::from_utf8(anf_bytes.clone()).map_err(|e| e.to_string())?;
        let manifest_bytes = read_bounded(manifest_path, 1_048_576)?;
        if sha256(&manifest_bytes) != expected_manifest {
            return Err("retained exporter manifest differs from frozen input pin".into());
        }
        let manifest: Value = serde_json::from_slice(&manifest_bytes).map_err(|e| e.to_string())?;
        validate_manifest(&manifest, point)?;
        let magma_bytes = read_bounded(&source_dir.join("instance.magma"), 2 * 1024 * 1024)?;
        for (role, name, bytes) in [
            ("wdsat_anf", "instance.anf", &anf_bytes),
            ("cryptominisat_xor_dimacs", "instance.xor.cnf", &cnf_bytes),
            ("magma_boolean_f4", "instance.magma", &magma_bytes),
        ] {
            let descriptor = &manifest["exports"][role];
            if descriptor["path"] != name
                || descriptor["bytes"] != bytes.len()
                || descriptor["blake3"] != blake3::hash(bytes).to_hex().as_str()
            {
                return Err(format!("retained {role} export bytes differ from manifest"));
            }
        }
        let width = manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
            .as_u64()
            .and_then(|v| usize::try_from(v).ok())
            .ok_or("missing source width")?;
        validate_source_syntax(&anf, &cnf, width)?;
        let preparation_bytes = read_bounded(preparation_path, 16 * 1024 * 1024)?;
        let preparation: Value =
            serde_json::from_slice(&preparation_bytes).map_err(|e| e.to_string())?;
        let state = PreparedState::load(&preparation)?;
        let public = FastPoint::affine(point[0], point[1]);
        if !state.curve.is_on_curve(public)
            || state.curve.mul_u64(public, ORDER) != FastPoint::INFINITY
        {
            return Err("disclosed public point fails independent curve/subgroup check".into());
        }
        let geometric = geometric_witness(&state, public);
        let geometric_class_valid = match expected_status {
            "SAT_MODEL" => geometric.is_some(),
            "SOURCE_UNSAT" => geometric.is_none(),
            _ => true,
        };
        let timeout = Duration::from_millis(timeout_ms);

        let cold_input = out.join("cold-input.xor.cnf");
        save_new(&cold_input, &cnf_bytes)?;
        let (mut cold, cold_argv) = launch(cms, Some(&cold_input), out, "cold", conflict_budget)?;
        let cold_pid = cold.0.id();
        let cold_result = wait(&mut cold, timeout)?;
        drop(cold);
        let cold = outcome(
            out,
            "cold",
            &cold_argv,
            cold_pid,
            cold_result,
            &anf,
            &cnf,
            &state,
            public,
        )?;

        let (mut prepared, prepared_argv) = launch(cms, None, out, "stdin", conflict_budget)?;
        let prepared_pid = prepared.0.id();
        let deadline = Instant::now() + timeout;
        let reader_ready_wait_ns = reader_ready(out, &mut prepared, deadline)?;
        let mut writer = prepared
            .0
            .stdin
            .take()
            .ok_or("missing prestarted CMS stdin")?;
        send_stdin(&mut writer, &cnf_bytes, &mut prepared, deadline)?;
        drop(writer);
        let remaining = deadline.saturating_duration_since(Instant::now());
        let prepared_result = wait(&mut prepared, remaining)?;
        drop(prepared);
        let stdin_run = outcome(
            out,
            "stdin",
            &prepared_argv,
            prepared_pid,
            prepared_result,
            &anf,
            &cnf,
            &state,
            public,
        )?;

        let cms_unchanged = sha256(&read_bounded(cms, 8 * 1024 * 1024)?) == cms_sha;
        let cold_input_unchanged =
            sha256(&read_bounded(&cold_input, MAX_CNF_BYTES)?) == sha256(&cnf_bytes);
        let source_unchanged = sha256(&read_bounded(cnf_path, MAX_CNF_BYTES)?)
            == sha256(&cnf_bytes)
            && sha256(&read_bounded(anf_path, MAX_ANF_BYTES)?) == sha256(&anf_bytes)
            && sha256(&read_bounded(manifest_path, 1_048_576)?) == sha256(&manifest_bytes)
            && sha256(&read_bounded(
                &source_dir.join("instance.magma"),
                2 * 1024 * 1024,
            )?) == sha256(&magma_bytes);
        let status_parity = cold["status"] == stdin_run["status"];
        let source_model_parity = cold["model_check"]["source_model_valid"]
            == stdin_run["model_check"]["source_model_valid"];
        Ok(
            json!({"schema_version":1,"question":"disclosed-cms-prestarted-stdin-control-v1",
        "cms_path":cms,"cnf_path":cnf_path,"anf_path":anf_path,
        "manifest_path":manifest_path,"preparation_path":preparation_path,"public_point":point,
        "preparation_sha256":sha256(&preparation_bytes),
        "geometric_witness":geometric,"geometric_class_valid":geometric_class_valid,
        "manifest_sha256":sha256(&manifest_bytes),
        "cms_sha256":cms_sha,"cms_unchanged":cms_unchanged,
        "cold_input_unchanged":cold_input_unchanged,"source_unchanged":source_unchanged,
        "cnf_sha256":sha256(&cnf_bytes),"anf_sha256":sha256(&anf_bytes),
        "magma_sha256":sha256(&magma_bytes),
        "conflict_budget":conflict_budget,"timeout_ms":timeout_ms,
        "expected_status":expected_status,
        "cold":cold,"stdin":stdin_run,"reader_ready_wait_ns":reader_ready_wait_ns,
        "reader_ready_before_cnf_delivery":true,
        "status_parity":status_parity,"source_model_parity":source_model_parity,
        "prestarted_stdin_compatible":cms_unchanged && cold_input_unchanged && source_unchanged && geometric_class_valid && status_parity
            && source_model_parity && cold["status"] == expected_status,
        "native_child_group_drain_audited":false,
        "source_bound_target_admitted":false,"online_speedup":null}),
        )
    }

    pub fn run_main() -> ExitCode {
        let argv = env::args_os().collect::<Vec<_>>();
        if argv.len() != 14 {
            eprintln!("usage: ic_sat_stdin_probe CMS CNF ANF MANIFEST PREPARATION POINT_X POINT_Y NEW_OUTPUT_DIR EXPECTED_CMS_SHA256 EXPECTED_MANIFEST_SHA256 CONFLICT_BUDGET TIMEOUT_MS EXPECTED_STATUS");
            return ExitCode::FAILURE;
        }
        let out = PathBuf::from(&argv[8]);
        if let Err(error) = fs::create_dir(&out) {
            eprintln!("new output directory required: {error}");
            return ExitCode::FAILURE;
        }
        let expected = argv[9].to_string_lossy();
        let expected_manifest = argv[10].to_string_lossy();
        let expected_status = argv[13].to_string_lossy();
        let point_x = argv[6].to_string_lossy().parse::<u64>();
        let point_y = argv[7].to_string_lossy().parse::<u64>();
        let conflict_budget = argv[11].to_string_lossy().parse::<u64>();
        let timeout_ms = argv[12].to_string_lossy().parse::<u64>();
        let result = point_x
            .map_err(|e| e.to_string())
            .and_then(|x| point_y.map_err(|e| e.to_string()).map(|y| [x, y]))
            .and_then(|point| {
                conflict_budget
                    .map_err(|e| e.to_string())
                    .map(|budget| (point, budget))
            })
            .and_then(|(point, budget)| {
                timeout_ms
                    .map_err(|e| e.to_string())
                    .map(|timeout| (point, budget, timeout))
            })
            .and_then(|(point, budget, timeout)| {
                run(
                    Path::new(&argv[1]),
                    Path::new(&argv[2]),
                    Path::new(&argv[3]),
                    Path::new(&argv[4]),
                    Path::new(&argv[5]),
                    point,
                    &out,
                    &expected,
                    &expected_manifest,
                    budget,
                    timeout,
                    &expected_status,
                )
            });
        let report = match &result {
            Ok(value) => value.clone(),
            Err(error) => {
                json!({"schema_version":1,"question":"disclosed-cms-prestarted-stdin-control-v1",
            "failure":error,"prestarted_stdin_compatible":false,"source_bound_target_admitted":false,
            "native_child_group_drain_audited":false,"online_speedup":null})
            }
        };
        let bytes = serde_json::to_vec_pretty(&report).expect("report serialisation");
        if let Err(error) = save_new(&out.join("result.json"), &bytes) {
            eprintln!("cannot retain probe result: {error}");
            return ExitCode::FAILURE;
        }
        println!("{}", out.join("result.json").display());
        if result
            .as_ref()
            .is_ok_and(|value| value["prestarted_stdin_compatible"] == true)
        {
            ExitCode::SUCCESS
        } else {
            ExitCode::FAILURE
        }
    }
}

#[cfg(unix)]
fn main() -> std::process::ExitCode {
    unix_probe::run_main()
}

#[cfg(not(unix))]
fn main() -> std::process::ExitCode {
    eprintln!("ic_sat_stdin_probe requires Unix child-process control");
    std::process::ExitCode::FAILURE
}

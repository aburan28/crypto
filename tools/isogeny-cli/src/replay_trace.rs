//! Opt-in, fresh per-replay JSONL evidence. No arithmetic or admission decisions.
use serde_json::json;
use std::{
    cell::RefCell,
    fs::{File, OpenOptions},
    io::{self, Write},
    path::Path,
    rc::Rc,
    time::Instant,
};
thread_local! { static ACTIVE: RefCell<Option<Rc<Trace>>> = const { RefCell::new(None) }; }
struct Frame {
    id: u64,
    parent: Option<u64>,
    name: String,
}
struct State {
    file: File,
    next_id: u64,
    stack: Vec<Frame>,
}
struct Trace {
    start: Instant,
    state: RefCell<State>,
}
impl Trace {
    fn record(
        &self,
        state: &mut State,
        frame: &Frame,
        event: &str,
        status: Option<&str>,
    ) -> io::Result<()> {
        serde_json::to_writer(
            &mut state.file,
            &json!({
                "schema":"isogeny-replay-stage/v1", "stage":frame.name, "event":event,
                "stage_id":frame.id, "parent_id":frame.parent, "pid":std::process::id(),
                "elapsed_ms_operational_only":self.start.elapsed().as_millis(), "status":status,
            }),
        )?;
        state.file.write_all(b"\n")?;
        state.file.flush()
    }
    fn begin(&self, name: String) -> io::Result<u64> {
        let mut state = self.state.borrow_mut();
        let frame = Frame {
            id: state.next_id,
            parent: state.stack.last().map(|f| f.id),
            name,
        };
        state.next_id += 1;
        self.record(&mut state, &frame, "begin", None)?;
        let id = frame.id;
        state.stack.push(frame);
        Ok(id)
    }
    fn end_id(&self, id: u64, event: &str, unwind: bool) -> io::Result<()> {
        let mut state = self.state.borrow_mut();
        loop {
            let frame = state
                .stack
                .pop()
                .expect("replay trace stage stack underflow");
            if frame.id == id {
                return self.record(&mut state, &frame, event, None);
            }
            assert!(unwind, "replay trace stage was not closed");
            self.record(&mut state, &frame, "error", None)?;
        }
    }
    fn end_name(&self, name: &str, event: &str) -> io::Result<()> {
        let id = {
            let state = self.state.borrow();
            let frame = state
                .stack
                .last()
                .expect("replay trace stage stack underflow");
            assert_eq!(frame.name, name, "replay trace stage ordering");
            frame.id
        };
        self.end_id(id, event, false)
    }
    fn terminal(&self, status: &str) -> io::Result<()> {
        let mut state = self.state.borrow_mut();
        self.record(
            &mut state,
            &Frame {
                id: 0,
                parent: None,
                name: "total".into(),
            },
            "terminal",
            Some(status),
        )
    }
}
fn current() -> Option<Rc<Trace>> {
    ACTIVE.with(|active| active.borrow().clone())
}
fn require(result: io::Result<()>) {
    result.expect("requested replay trace could not be written/flushed");
}
/// A total session includes dispatcher input parsing and the selected field engine.
pub struct Session {
    trace: Option<Rc<Trace>>,
    total_id: Option<u64>,
    root: bool,
}
impl Session {
    pub fn start() -> Self {
        if let Some(trace) = current() {
            return Self {
                trace: Some(trace),
                total_id: None,
                root: false,
            };
        }
        let Some(path) = std::env::var_os("ISOGENY_REPLAY_TRACE") else {
            return Self {
                trace: None,
                total_id: None,
                root: false,
            };
        };
        let path = Path::new(&path);
        assert!(
            path.is_absolute(),
            "requested replay trace needs an absolute fresh path"
        );
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .expect(
                "requested replay trace needs an absent file and an available parent directory",
            );
        let trace = Rc::new(Trace {
            start: Instant::now(),
            state: RefCell::new(State {
                file,
                next_id: 0,
                stack: vec![],
            }),
        });
        let id = trace
            .begin("total".into())
            .expect("requested replay trace could not begin");
        ACTIVE.with(|active| *active.borrow_mut() = Some(trace.clone()));
        Self {
            trace: Some(trace),
            total_id: Some(id),
            root: true,
        }
    }
    pub fn finish(mut self) {
        if self.root {
            let trace = self.trace.as_ref().unwrap();
            let id = self.total_id.take().unwrap();
            require(trace.end_id(id, "end", false));
            require(trace.terminal("COMPLETED"));
            ACTIVE.with(|active| *active.borrow_mut() = None);
            self.root = false;
        }
    }
}
impl Drop for Session {
    fn drop(&mut self) {
        if self.root {
            let trace = self.trace.as_ref().unwrap();
            if let Some(id) = self.total_id {
                if let Err(error) = trace.end_id(id, "error", true) {
                    eprintln!("requested replay trace failed during error recording: {error}");
                }
            }
            if let Err(error) = trace.terminal("ERROR") {
                eprintln!("requested replay trace failed during terminal recording: {error}");
            }
            ACTIVE.with(|active| *active.borrow_mut() = None);
        }
    }
}
pub struct Stage {
    trace: Option<Rc<Trace>>,
    id: Option<u64>,
}
pub fn stage(name: impl Into<String>) -> Stage {
    let trace = current();
    let id = trace.as_ref().map(|trace| {
        trace
            .begin(name.into())
            .expect("requested replay trace could not begin stage")
    });
    Stage { trace, id }
}
impl Stage {
    pub fn finish(mut self) {
        if let Some(id) = self.id.take() {
            require(self.trace.as_ref().unwrap().end_id(id, "end", false));
            self.id = None;
        }
    }
}
impl Drop for Stage {
    fn drop(&mut self) {
        if let Some(id) = self.id {
            if let Err(error) = self.trace.as_ref().unwrap().end_id(id, "error", true) {
                eprintln!("requested replay trace failed during stage error recording: {error}");
            }
        }
    }
}
/// Observer bridge used by the independent kernel verifier's no-op-compatible API.
pub fn event(name: &'static str, event: &'static str) {
    if let Some(trace) = current() {
        match event {
            "begin" => {
                trace
                    .begin(name.into())
                    .expect("requested replay trace could not begin kernel stage");
            }
            "end" | "error" => require(trace.end_name(name, event)),
            _ => panic!("unknown internal replay trace event"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs,
        process::{Command, Stdio},
        thread,
        time::{Duration, SystemTime},
    };
    #[test]
    #[ignore = "child fixture selected by the interruption control"]
    fn interruption_worker() {
        if std::env::var_os("ISOGENY_TRACE_CONTROL_CHILD").is_none() {
            return;
        }
        let session = Session::start();
        let pending = stage("interruption_fixture_wait");
        thread::sleep(Duration::from_secs(10));
        pending.finish();
        session.finish();
    }
    #[test]
    fn flushed_begin_survives_bounded_child_interruption() {
        let retained = std::env::var_os("ISOGENY_TRACE_UNIT_EVIDENCE");
        let root = retained
            .as_ref()
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|| {
                std::env::temp_dir().join(format!(
                    "isogeny-trace-interruption-{}-{}",
                    std::process::id(),
                    SystemTime::now()
                        .duration_since(SystemTime::UNIX_EPOCH)
                        .unwrap()
                        .as_nanos()
                ))
            });
        fs::create_dir(&root).unwrap();
        let root = fs::canonicalize(root).unwrap();
        let path = root.join("interrupted.jsonl");
        let mut child = Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "replay_trace::tests::interruption_worker",
                "--ignored",
                "--nocapture",
            ])
            .env("ISOGENY_TRACE_CONTROL_CHILD", "1")
            .env("ISOGENY_REPLAY_TRACE", &path)
            .stdout(Stdio::from(
                File::create(root.join("child.stdout.log")).unwrap(),
            ))
            .stderr(Stdio::from(
                File::create(root.join("child.stderr.log")).unwrap(),
            ))
            .spawn()
            .unwrap();
        let started = Instant::now();
        let mut observed = false;
        while started.elapsed() < Duration::from_secs(5) {
            if let Ok(text) = fs::read_to_string(&path) {
                observed = text
                    .lines()
                    .filter_map(|line| serde_json::from_str::<serde_json::Value>(line).ok())
                    .any(|record| {
                        record["stage"] == "interruption_fixture_wait" && record["event"] == "begin"
                    });
                if observed {
                    break;
                }
            }
            thread::sleep(Duration::from_millis(5));
        }
        let _ = child.kill();
        let status = child.wait().unwrap();
        assert!(
            observed,
            "flushed begin must be visible before interrupting"
        );
        let records: Vec<serde_json::Value> = fs::read_to_string(&path)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        assert_eq!(records.len(), 2);
        assert!(records.iter().all(|record| record["event"] == "begin"));
        assert_eq!(records[0]["stage"], "total");
        assert_eq!(records[1]["stage"], "interruption_fixture_wait");
        fs::write(root.join("control.json"), serde_json::to_vec_pretty(&json!({"schema":"isogeny-trace-interruption-control/v1", "cryptographic_work":false, "flushed_begin_observed_before_kill":observed, "child_exit_code":status.code(), "unfinished_stage_marked_complete":false})).unwrap()).unwrap();
        if retained.is_none() {
            fs::remove_dir_all(root).unwrap();
        }
    }
}

#[cfg(test)]
mod write_failure_tests {
    use super::*;
    #[test]
    fn requested_trace_write_failure_propagates_without_double_panic() {
        let path = std::env::temp_dir().join(format!(
            "isogeny-trace-readonly-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, b"retained input\n").unwrap();
        // Read-only descriptor deterministically rejects writes, including under
        // privileged test users. No production telemetry/trace override exists.
        let trace = Rc::new(Trace {
            start: Instant::now(),
            state: RefCell::new(State {
                file: File::open(&path).unwrap(),
                next_id: 1,
                stack: vec![Frame {
                    id: 0,
                    parent: None,
                    name: "total".into(),
                }],
            }),
        });
        ACTIVE.with(|active| *active.borrow_mut() = Some(trace.clone()));
        let session = Session {
            trace: Some(trace),
            total_id: Some(0),
            root: true,
        };
        let error = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| session.finish()));
        assert!(error.is_err());
        assert!(current().is_none());
        assert_eq!(std::fs::read(&path).unwrap(), b"retained input\n");
        std::fs::remove_file(path).unwrap();
    }
}

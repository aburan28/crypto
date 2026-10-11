//! `spotlab demo-job`: a resumable reference workload.
//!
//! It extends a SHA-256 hash chain `h_i = SHA256(h_{i-1} || i)` from `h_0 =
//! SHA256("spotlab")` to `--to N`, saving `(i, h_i)` to
//! `$SPOTLAB_CHECKPOINT_DIR/state.json` every `--every` steps and on SIGTERM,
//! and restoring from it at start-up. Because each step depends on the last,
//! a final hash that matches the uninterrupted chain shows that every resume
//! continued from exactly where the previous attempt stopped.

use anyhow::Result;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::util::write_atomic;

static TERM: AtomicBool = AtomicBool::new(false);

extern "C" fn on_term(_: libc::c_int) {
    TERM.store(true, Ordering::SeqCst);
}

#[derive(Serialize, Deserialize)]
struct State {
    i: u64,
    h: String,
    starts: u32,
}

pub fn step(h: &[u8; 32], i: u64) -> [u8; 32] {
    let mut d = Sha256::new();
    d.update(h);
    d.update(i.to_le_bytes());
    d.finalize().into()
}

#[cfg(test)]
pub fn chain(n: u64) -> [u8; 32] {
    let mut h: [u8; 32] = Sha256::digest(b"spotlab").into();
    for i in 1..=n {
        h = step(&h, i);
    }
    h
}

pub fn run(to: u64, every: u64, step_us: u64) -> Result<()> {
    unsafe {
        libc::signal(libc::SIGTERM, on_term as *const () as libc::sighandler_t);
    }
    let ck = PathBuf::from(
        std::env::var("SPOTLAB_CHECKPOINT_DIR").unwrap_or_else(|_| "checkpoint".into()),
    );
    let out =
        PathBuf::from(std::env::var("SPOTLAB_OUTPUT_DIR").unwrap_or_else(|_| "output".into()));
    let state_path = ck.join("state.json");
    let mut st = match std::fs::read(&state_path) {
        Ok(b) => serde_json::from_slice::<State>(&b)?,
        Err(_) => State {
            i: 0,
            h: hex::encode(Sha256::digest(b"spotlab")),
            starts: 0,
        },
    };
    st.starts += 1;
    eprintln!("demo-job: start {} at step {}", st.starts, st.i);
    let mut h: [u8; 32] = hex::decode(&st.h)?
        .try_into()
        .map_err(|_| anyhow::anyhow!("bad hash in state"))?;
    let save = |i: u64, h: &[u8; 32], starts: u32| -> Result<()> {
        let s = State {
            i,
            h: hex::encode(h),
            starts,
        };
        write_atomic(&state_path, &serde_json::to_vec(&s)?)
    };
    while st.i < to {
        if TERM.load(Ordering::SeqCst) {
            save(st.i, &h, st.starts)?;
            eprintln!("demo-job: SIGTERM at step {}; saved", st.i);
            std::process::exit(143);
        }
        st.i += 1;
        h = step(&h, st.i);
        if step_us > 0 {
            std::thread::sleep(std::time::Duration::from_micros(step_us));
        }
        if st.i % every.max(1) == 0 {
            save(st.i, &h, st.starts)?;
        }
    }
    save(st.i, &h, st.starts)?;
    let result = serde_json::json!({"steps": st.i, "hash": hex::encode(h), "starts": st.starts});
    write_atomic(
        &out.join("result.json"),
        serde_json::to_string_pretty(&result)?.as_bytes(),
    )?;
    println!("{result}");
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn chain_is_deterministic_and_stepwise() {
        let a = super::chain(10);
        assert_eq!(super::step(&super::chain(9), 10), a);
        assert_ne!(super::chain(11), a);
    }
}

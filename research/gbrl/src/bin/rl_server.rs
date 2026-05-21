// JSON-RPC over stdio. Python trainer spawns this binary, sends one request
// per line, reads one response per line.
//
// Request shapes:
//   {"method": "Reset",   "instance": "cyclic-5"}
//   {"method": "Step",    "pair_idx": 3}
//   {"method": "Observe"}
//
// Response shape:
//   {"ok": bool, "obs": <Observation> | null, "reward": f64 | null,
//    "done": bool | null, "error": string | null}
//
//     cargo run --release --bin rl_server

use gbrl::env::Env;
use gbrl::semaev::{cyclic, semaev_decomposition_3_instance, semaev_decomposition_diverse,
                    semaev_decomposition_instance};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::io::{self, BufRead, Write};

#[derive(Deserialize)]
#[serde(tag = "method")]
enum Request {
    Reset { instance: String },
    Step { pair_idx: usize },
    Observe,
}

#[derive(Serialize)]
struct Response {
    ok: bool,
    obs: Option<Value>,
    reward: Option<f64>,
    done: Option<bool>,
    error: Option<String>,
}

impl Response {
    fn err(msg: impl Into<String>) -> Self {
        Response { ok: false, obs: None, reward: None, done: None, error: Some(msg.into()) }
    }
}

fn make_instance(name: &str) -> Option<Env> {
    if let Some(rest) = name.strip_prefix("cyclic-") {
        let n: usize = rest.parse().ok()?;
        return Some(Env::new(cyclic(n), 200_000));
    }
    // Semaev 2-decomposition system. The number is the factor-base size,
    // which controls the leading-monomial structure (degree of F polys).
    // e.g. "semaev2-8" is the same family used in decomp_demo.
    if let Some(rest) = name.strip_prefix("semaev2d-") {
        let seed: u64 = rest.parse().ok()?;
        return Some(Env::new(semaev_decomposition_diverse(seed), 200_000));
    }
    if let Some(rest) = name.strip_prefix("semaev2-") {
        let n: usize = rest.parse().ok()?;
        return Some(Env::new(semaev_decomposition_instance(n), 200_000));
    }
    // 3-decomposition system (much harder than m=2). Number is factor-base size.
    if let Some(rest) = name.strip_prefix("semaev3-") {
        let n: usize = rest.parse().ok()?;
        return Some(Env::new(semaev_decomposition_3_instance(n), 200_000));
    }
    None
}

fn handle(req: Request, env: &mut Option<Env>) -> Response {
    match req {
        Request::Reset { instance } => match make_instance(&instance) {
            None => Response::err(format!("unknown instance: {}", instance)),
            Some(e) => {
                *env = Some(e);
                let obs = env.as_ref().unwrap().observe();
                let done = obs.done;
                Response {
                    ok: true,
                    obs: Some(serde_json::to_value(obs).unwrap()),
                    reward: None,
                    done: Some(done),
                    error: None,
                }
            }
        },
        Request::Step { pair_idx } => match env.as_mut() {
            None => Response::err("no env; call Reset first"),
            Some(e) => {
                if pair_idx >= e.state.pairs.len() {
                    return Response::err(format!(
                        "pair_idx {} out of range (have {} pairs)",
                        pair_idx, e.state.pairs.len()
                    ));
                }
                let (r, d) = e.step(pair_idx);
                let obs = e.observe();
                Response {
                    ok: true,
                    obs: Some(serde_json::to_value(obs).unwrap()),
                    reward: Some(r),
                    done: Some(d),
                    error: None,
                }
            }
        },
        Request::Observe => match env.as_ref() {
            None => Response::err("no env; call Reset first"),
            Some(e) => {
                let obs = e.observe();
                let done = obs.done;
                Response {
                    ok: true,
                    obs: Some(serde_json::to_value(obs).unwrap()),
                    reward: None,
                    done: Some(done),
                    error: None,
                }
            }
        },
    }
}

fn main() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut env: Option<Env> = None;

    for line in stdin.lock().lines() {
        let line = match line { Ok(l) => l, Err(_) => break };
        if line.trim().is_empty() { continue; }

        let resp = match serde_json::from_str::<Request>(&line) {
            Err(e) => Response::err(format!("parse: {}", e)),
            Ok(req) => handle(req, &mut env),
        };

        let mut h = stdout.lock();
        writeln!(h, "{}", serde_json::to_string(&resp).unwrap()).unwrap();
        h.flush().unwrap();
    }
}

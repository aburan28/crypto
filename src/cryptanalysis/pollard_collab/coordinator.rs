//! Coordinator transport — one hub with a URL, agents that dial out.
//!
//! The gossip transport in [`net`](super::net) assumes every node can
//! *accept* a connection.  On a cloud fleet that is rarely true: the
//! agents sit behind NAT, in private subnets, on laptops, or in
//! containers with no inbound rule, while exactly one host — an EC2
//! instance with an elastic IP or a load balancer in front of it — is
//! reachable by everybody.  This module is that shape:
//!
//! ```text
//!   agent (private subnet)  ──────┐
//!   agent (laptop, NAT)     ──────┼──►  https://coordinator:8080   (EC2)
//!   agent (spot fleet)      ──────┘        GET /v1/channel
//!                                          Upgrade: rho-collab/1
//! ```
//!
//! Every connection is **opened by the agent**.  Once it is up the
//! coordinator uses that same socket in the other direction — the
//! *reverse channel* — to push check-ins, other agents' distinguished
//! points and the solution down to an agent it could never have
//! dialled.  Nothing is polled: an agent learns about a collision
//! found on another continent as soon as the hub has merged it.
//!
//! ## What the coordinator is, and is not
//!
//! It is a *rendezvous*, not an authority.  It holds the same
//! [`SharedState`] CRDT every agent holds, merges what arrives with
//! the same verification (`apply(.., verify = true)`), and hands out
//! nothing: agents still choose their own units by the leasing rules
//! in [`state`](super::state).  Losing the coordinator loses no work —
//! agents keep walking, keep their own DP tables, and reconverge when
//! it comes back, exactly as after a network partition.  What the hub
//! adds over a mesh is reachability and one place to look.
//!
//! ## The wire
//!
//! HTTP/1.1, so it can sit behind an ALB, nginx or Caddy for TLS:
//!
//! | route | method | purpose |
//! |---|---|---|
//! | `/healthz` | GET | load-balancer health check |
//! | `/v1/job` | GET | the job document, so an agent needs only the URL |
//! | `/v1/status` | GET | [`Progress`] as JSON |
//! | `/v1/sync` | POST | one-shot pull+push, for `status` and for agents that cannot hold a socket open |
//! | `/v1/channel` | GET + `Upgrade: rho-collab/1` | the reverse channel |
//!
//! The channel handshake is an ordinary HTTP upgrade (`101 Switching
//! Protocols`), the same move WebSocket makes, which is what lets it
//! through proxies.  After it, both directions speak one JSON
//! [`Frame`] per line:
//!
//! ```text
//! agent → hub   hello {job_id, peer, known}      who I am, what I have
//! hub   → agent batch {checkins, known, solution}  what you lack
//! agent → hub   push  {checkins}                   what I found
//! hub   → agent ack   {accepted, rejected}
//! hub   → agent batch {…}                          unprompted, whenever others report
//! agent → hub   pull  {known}                      resync after a gap
//! either        ping  {}                           keepalive through idle-timeout proxies
//! ```
//!
//! ## Authentication
//!
//! A coordinator on a public address takes a bearer token
//! (`Authorization: Bearer …`), compared in constant time.  The token
//! gates *writes as much as reads*: an unauthenticated peer could
//! otherwise flood the log with (verifiable but useless) check-ins.
//! It is a shared secret, not an identity — the trust model of
//! `POLLARD_COLLAB_DESIGN.md` §10 is unchanged, and every DP is still
//! verified on arrival, so the token buys access control, not
//! integrity.  Plain HTTP carries it, so terminate TLS in front of the
//! process; `deploy/rho-coordinator/` does exactly that.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream, ToSocketAddrs};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use serde::{Deserialize, Serialize};
use subtle::ConstantTimeEq;

use super::job::{hex_of, JobContext, JobSpec};
use super::state::{now_secs, CheckIn, Progress, SharedState};

/// Protocol token in the `Upgrade:` header.
pub const CHANNEL_PROTOCOL: &str = "rho-collab/1";

/// Cap on one frame or request body (a batch of check-ins).
const MAX_BODY: usize = 64 << 20;

/// Environment variables the CLI falls back to, so a fleet image can
/// be baked once and pointed at a hub by the instance's environment.
pub const URL_ENV: &str = "RHO_COORDINATOR_URL";
/// Environment variable holding the bearer token.
pub const TOKEN_ENV: &str = "RHO_COORDINATOR_TOKEN";

/// One line on the reverse channel.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Frame {
    /// First frame from the agent: who it is and what it already has.
    Hello {
        job_id: String,
        peer: String,
        known: BTreeMap<String, u64>,
    },
    /// "Send me what I lack" — also re-arms the hub's push cursor.
    Pull { known: BTreeMap<String, u64> },
    /// Check-ins the receiver was missing, plus the sender's vector.
    Batch {
        checkins: Vec<CheckIn>,
        known: BTreeMap<String, u64>,
        solution: Option<String>,
    },
    /// Check-ins the sender believes the receiver lacks.
    Push { checkins: Vec<CheckIn> },
    /// Result of applying a [`Frame::Push`].
    Ack { accepted: usize, rejected: usize },
    /// Keepalive.  Idle NAT and proxy timeouts are the reason the
    /// reverse channel needs one at all.
    Ping,
    /// Fatal for the connection; the sender hangs up after it.
    Error { message: String },
}

fn write_frame(w: &mut impl Write, f: &Frame) -> std::io::Result<()> {
    let mut line = serde_json::to_vec(f)?;
    line.push(b'\n');
    w.write_all(&line)?;
    w.flush()
}

fn read_frame(r: &mut impl BufRead) -> std::io::Result<Option<Frame>> {
    let mut line = String::new();
    if r.read_line(&mut line)? == 0 {
        return Ok(None);
    }
    if line.len() > MAX_BODY {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "frame too long",
        ));
    }
    if line.trim().is_empty() {
        return Ok(Some(Frame::Ping));
    }
    serde_json::from_str(&line)
        .map(Some)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

fn token_ok(expected: &Option<String>, presented: Option<&str>) -> bool {
    match expected {
        None => true,
        Some(t) => match presented {
            None => false,
            Some(p) => p.as_bytes().ct_eq(t.as_bytes()).into(),
        },
    }
}

// ── server ──────────────────────────────────────────────────────────────────

/// Coordinator tuning.
#[derive(Clone, Debug)]
pub struct CoordinatorConfig {
    /// Bearer token agents must present; `None` disables auth (use
    /// only on a private subnet or a loopback bind).
    pub token: Option<String>,
    /// How often an idle channel is examined for check-ins to push,
    /// and how often a keepalive goes out.
    pub push_interval: Duration,
    /// Seconds of silence after which a channel is dropped.  A lane
    /// that is merely slow keeps its units: this closes sockets, not
    /// leases.
    pub idle_timeout: Duration,
}

impl Default for CoordinatorConfig {
    fn default() -> Self {
        Self {
            token: None,
            push_interval: Duration::from_millis(500),
            idle_timeout: Duration::from_secs(300),
        }
    }
}

/// Live counters, for `/v1/status` and the operator's log line.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct HubStats {
    /// Channels open right now.
    pub agents: u64,
    /// Channels opened since start.
    pub channels_total: u64,
    /// Check-ins accepted from agents.
    pub accepted: u64,
    /// Check-ins or DPs refused (bad job, forged point, bad token).
    pub rejected: u64,
    /// Check-ins pushed down reverse channels.
    pub pushed: u64,
    /// Requests refused for a bad or missing token.
    pub unauthorized: u64,
}

#[derive(Default)]
struct Counters {
    agents: AtomicU64,
    channels_total: AtomicU64,
    accepted: AtomicU64,
    rejected: AtomicU64,
    pushed: AtomicU64,
    unauthorized: AtomicU64,
}

impl Counters {
    fn snapshot(&self) -> HubStats {
        HubStats {
            agents: self.agents.load(Ordering::Relaxed),
            channels_total: self.channels_total.load(Ordering::Relaxed),
            accepted: self.accepted.load(Ordering::Relaxed),
            rejected: self.rejected.load(Ordering::Relaxed),
            pushed: self.pushed.load(Ordering::Relaxed),
            unauthorized: self.unauthorized.load(Ordering::Relaxed),
        }
    }
}

/// Called with every check-in the hub accepts, after it is merged —
/// the durability hook (write it to a mailbox directory, ship it to
/// S3, …).  See [`Coordinator::start`].
pub type OnCheckIn = Box<dyn Fn(&CheckIn) + Send + Sync>;

struct Shared {
    ctx: Arc<JobContext>,
    spec_json: String,
    state: Arc<Mutex<SharedState>>,
    cfg: CoordinatorConfig,
    counters: Counters,
    on_checkin: Option<OnCheckIn>,
    lease_secs: u64,
}

impl Shared {
    /// Merge one check-in, count it, and run the durability hook.
    fn absorb(&self, ci: &CheckIn) -> (usize, usize) {
        let outcome = {
            let mut st = self.state.lock().unwrap();
            st.apply(&self.ctx, ci, now_secs(), true)
        };
        match outcome {
            Ok(o) => {
                if o.new {
                    self.counters.accepted.fetch_add(1, Ordering::Relaxed);
                    if let Some(hook) = &self.on_checkin {
                        hook(ci);
                    }
                }
                self.counters
                    .rejected
                    .fetch_add(o.rejected_dps as u64, Ordering::Relaxed);
                (o.accepted_dps, o.rejected_dps)
            }
            Err(_) => {
                self.counters.rejected.fetch_add(1, Ordering::Relaxed);
                (0, 1)
            }
        }
    }
}

/// A running coordinator.  Drop it or call [`stop`](Self::stop).
pub struct Coordinator {
    addr: std::net::SocketAddr,
    shared: Arc<Shared>,
    stop: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
}

impl Coordinator {
    /// Bind `addr` (port 0 for an ephemeral one) and serve in the
    /// background.  `on_checkin` runs for every check-in newly
    /// accepted from an agent; use it to mirror the log somewhere that
    /// outlives the instance.
    pub fn start(
        addr: impl ToSocketAddrs,
        ctx: Arc<JobContext>,
        spec: &JobSpec,
        state: Arc<Mutex<SharedState>>,
        cfg: CoordinatorConfig,
        lease_secs: u64,
        on_checkin: Option<OnCheckIn>,
    ) -> std::io::Result<Self> {
        let listener = TcpListener::bind(addr)?;
        listener.set_nonblocking(true)?;
        let local = listener.local_addr()?;
        let shared = Arc::new(Shared {
            ctx,
            spec_json: spec.to_json(),
            state,
            cfg,
            counters: Counters::default(),
            on_checkin,
            lease_secs,
        });
        let stop = Arc::new(AtomicBool::new(false));
        let (s2, sh2) = (Arc::clone(&stop), Arc::clone(&shared));
        let handle = std::thread::spawn(move || {
            while !s2.load(Ordering::Relaxed) {
                match listener.accept() {
                    Ok((stream, _)) => {
                        let _ = stream.set_nonblocking(false);
                        let sh = Arc::clone(&sh2);
                        let st = Arc::clone(&s2);
                        std::thread::spawn(move || {
                            let _ = serve_connection(stream, &sh, &st);
                        });
                    }
                    Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                        std::thread::sleep(Duration::from_millis(20));
                    }
                    Err(_) => break,
                }
            }
        });
        Ok(Self {
            addr: local,
            shared,
            stop,
            handle: Some(handle),
        })
    }

    /// The bound address (useful when port 0 was requested).
    pub fn local_addr(&self) -> std::net::SocketAddr {
        self.addr
    }

    /// `http://<addr>`, the value agents pass as `--coordinator`.
    pub fn url(&self) -> String {
        format!("http://{}", self.addr)
    }

    pub fn stats(&self) -> HubStats {
        self.shared.counters.snapshot()
    }

    pub fn stop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

impl Drop for Coordinator {
    fn drop(&mut self) {
        self.stop();
    }
}

struct Request {
    method: String,
    path: String,
    authorization: Option<String>,
    upgrade: Option<String>,
    body: Vec<u8>,
}

fn read_request(reader: &mut BufReader<TcpStream>) -> std::io::Result<Option<Request>> {
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Ok(None);
    }
    let mut parts = line.split_whitespace();
    let method = parts.next().unwrap_or_default().to_string();
    let path = parts.next().unwrap_or_default().to_string();
    let mut authorization = None;
    let mut upgrade = None;
    let mut length = 0usize;
    loop {
        let mut h = String::new();
        let n = reader.read_line(&mut h)?;
        if n == 0 || h == "\r\n" || h == "\n" {
            break;
        }
        if let Some((name, value)) = h.split_once(':') {
            let name = name.trim().to_ascii_lowercase();
            let value = value.trim().to_string();
            match name.as_str() {
                "authorization" => authorization = Some(value),
                "upgrade" => upgrade = Some(value),
                "content-length" => length = value.parse().unwrap_or(0),
                _ => {}
            }
        }
    }
    if length > MAX_BODY {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "body too long",
        ));
    }
    let mut body = vec![0u8; length];
    if length > 0 {
        reader.read_exact(&mut body)?;
    }
    Ok(Some(Request {
        method,
        path,
        authorization,
        upgrade,
        body,
    }))
}

fn respond(w: &mut impl Write, status: u16, reason: &str, body: &str) -> std::io::Result<()> {
    let head = format!(
        "HTTP/1.1 {status} {reason}\r\nContent-Type: application/json\r\n\
         Content-Length: {}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    w.write_all(head.as_bytes())?;
    w.write_all(body.as_bytes())?;
    w.flush()
}

fn bearer(req: &Request) -> Option<&str> {
    req.authorization
        .as_deref()
        .and_then(|v| v.strip_prefix("Bearer "))
        .map(str::trim)
}

fn serve_connection(
    stream: TcpStream,
    shared: &Arc<Shared>,
    stop: &Arc<AtomicBool>,
) -> std::io::Result<()> {
    stream.set_read_timeout(Some(Duration::from_secs(30)))?;
    let mut writer = stream.try_clone()?;
    let mut reader = BufReader::new(stream);
    let req = match read_request(&mut reader)? {
        Some(r) => r,
        None => return Ok(()),
    };
    // `/healthz` is deliberately outside the token: a load balancer
    // holds no credential, and the reply says nothing about the job.
    if req.path.starts_with("/healthz") {
        return respond(&mut writer, 200, "OK", "{\"ok\":true}");
    }
    if !token_ok(&shared.cfg.token, bearer(&req)) {
        shared.counters.unauthorized.fetch_add(1, Ordering::Relaxed);
        return respond(
            &mut writer,
            401,
            "Unauthorized",
            "{\"error\":\"bad or missing bearer token\"}",
        );
    }
    match (
        req.method.as_str(),
        req.path.split('?').next().unwrap_or(""),
    ) {
        ("GET", "/v1/job") => respond(&mut writer, 200, "OK", &shared.spec_json),
        ("GET", "/v1/status") => {
            let body = serde_json::to_string_pretty(&hub_status(shared)).unwrap_or_default();
            respond(&mut writer, 200, "OK", &body)
        }
        ("POST", "/v1/sync") => {
            let frame: Frame = match serde_json::from_slice(&req.body) {
                Ok(f) => f,
                Err(e) => {
                    return respond(
                        &mut writer,
                        400,
                        "Bad Request",
                        &serde_json::json!({ "error": e.to_string() }).to_string(),
                    )
                }
            };
            let reply = one_shot(shared, frame);
            let body = serde_json::to_string(&reply).unwrap_or_default();
            respond(&mut writer, 200, "OK", &body)
        }
        ("GET", "/v1/channel") => {
            let wants = req
                .upgrade
                .as_deref()
                .map(|u| u.eq_ignore_ascii_case(CHANNEL_PROTOCOL))
                .unwrap_or(false);
            if !wants {
                return respond(
                    &mut writer,
                    426,
                    "Upgrade Required",
                    &serde_json::json!({ "error": format!("upgrade to {CHANNEL_PROTOCOL}") })
                        .to_string(),
                );
            }
            writer.write_all(
                format!(
                    "HTTP/1.1 101 Switching Protocols\r\nUpgrade: {CHANNEL_PROTOCOL}\r\n\
                     Connection: Upgrade\r\n\r\n"
                )
                .as_bytes(),
            )?;
            writer.flush()?;
            run_channel(reader, writer, shared, stop)
        }
        _ => respond(
            &mut writer,
            404,
            "Not Found",
            "{\"error\":\"no such route\"}",
        ),
    }
}

/// Everything `/v1/status` reports: job progress plus hub counters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HubStatus {
    pub job_id: String,
    pub progress: Progress,
    pub hub: HubStats,
}

fn hub_status(shared: &Arc<Shared>) -> HubStatus {
    let st = shared.state.lock().unwrap();
    HubStatus {
        job_id: shared.ctx.job_id.clone(),
        progress: st.progress(&shared.ctx, now_secs(), shared.lease_secs),
        hub: shared.counters.snapshot(),
    }
}

/// `/v1/sync`: apply what came in, answer with what the caller lacks.
fn one_shot(shared: &Arc<Shared>, frame: Frame) -> Frame {
    let (known, incoming) = match frame {
        Frame::Hello { job_id, known, .. } => {
            if job_id != shared.ctx.job_id {
                return Frame::Error {
                    message: "different job".into(),
                };
            }
            (known, Vec::new())
        }
        Frame::Pull { known } => (known, Vec::new()),
        Frame::Push { checkins } => (BTreeMap::new(), checkins),
        _ => {
            return Frame::Error {
                message: "expected hello, pull or push".into(),
            }
        }
    };
    for ci in &incoming {
        shared.absorb(ci);
    }
    let st = shared.state.lock().unwrap();
    Frame::Batch {
        checkins: st.delta_for(&known),
        known: st.version_vector(),
        solution: st.solution.as_ref().map(hex_of),
    }
}

/// The reverse channel, hub side.
///
/// Two threads: this one reads the agent's frames, a spawned one
/// pushes down whatever the agent is missing.  They share the socket's
/// write half behind a mutex and the agent's version vector behind
/// another, which the reader advances from `hello`/`pull` and the
/// writer advances by what it has sent.
fn run_channel(
    mut reader: BufReader<TcpStream>,
    writer: TcpStream,
    shared: &Arc<Shared>,
    stop: &Arc<AtomicBool>,
) -> std::io::Result<()> {
    // Blocking reads must not kill an agent that is merely quiet, so
    // the read timeout is a tick and the idle check is explicit.
    reader
        .get_ref()
        .set_read_timeout(Some(Duration::from_secs(5)))?;
    let out = Arc::new(Mutex::new(writer));
    let known = Arc::new(Mutex::new(BTreeMap::<String, u64>::new()));
    let alive = Arc::new(AtomicBool::new(true));
    shared.counters.agents.fetch_add(1, Ordering::Relaxed);
    shared
        .counters
        .channels_total
        .fetch_add(1, Ordering::Relaxed);

    let pusher = {
        let (shared, out, known, alive, stop) = (
            Arc::clone(shared),
            Arc::clone(&out),
            Arc::clone(&known),
            Arc::clone(&alive),
            Arc::clone(stop),
        );
        std::thread::spawn(move || {
            let mut last_ping = Instant::now();
            while alive.load(Ordering::Relaxed) && !stop.load(Ordering::Relaxed) {
                std::thread::sleep(shared.cfg.push_interval);
                let (delta, vector, solution) = {
                    let have = known.lock().unwrap().clone();
                    let st = shared.state.lock().unwrap();
                    (
                        st.delta_for(&have),
                        st.version_vector(),
                        st.solution.as_ref().map(hex_of),
                    )
                };
                if delta.is_empty() {
                    if last_ping.elapsed() >= shared.cfg.idle_timeout / 3 {
                        last_ping = Instant::now();
                        if write_frame(&mut *out.lock().unwrap(), &Frame::Ping).is_err() {
                            alive.store(false, Ordering::Relaxed);
                        }
                    }
                    continue;
                }
                let n = delta.len() as u64;
                let frame = Frame::Batch {
                    checkins: delta,
                    known: vector.clone(),
                    solution,
                };
                if write_frame(&mut *out.lock().unwrap(), &frame).is_err() {
                    alive.store(false, Ordering::Relaxed);
                    break;
                }
                shared.counters.pushed.fetch_add(n, Ordering::Relaxed);
                last_ping = Instant::now();
                // The agent now has at least what we just sent; a
                // `pull` from it can still lower this, which only
                // costs a resend.
                *known.lock().unwrap() = vector;
            }
        })
    };

    let mut last_seen = Instant::now();
    let result = loop {
        if stop.load(Ordering::Relaxed) || !alive.load(Ordering::Relaxed) {
            break Ok(());
        }
        match read_frame(&mut reader) {
            Ok(None) => break Ok(()),
            Ok(Some(frame)) => {
                last_seen = Instant::now();
                match frame {
                    Frame::Hello {
                        job_id,
                        known: theirs,
                        ..
                    } => {
                        if job_id != shared.ctx.job_id {
                            let _ = write_frame(
                                &mut *out.lock().unwrap(),
                                &Frame::Error {
                                    message: "different job".into(),
                                },
                            );
                            break Ok(());
                        }
                        *known.lock().unwrap() = theirs;
                    }
                    Frame::Pull { known: theirs } => *known.lock().unwrap() = theirs,
                    Frame::Push { checkins } => {
                        let (mut acc, mut rej) = (0, 0);
                        for ci in &checkins {
                            let (a, r) = shared.absorb(ci);
                            acc += a;
                            rej += r;
                        }
                        let _ = write_frame(
                            &mut *out.lock().unwrap(),
                            &Frame::Ack {
                                accepted: acc,
                                rejected: rej,
                            },
                        );
                    }
                    Frame::Ping => {}
                    _ => {}
                }
            }
            Err(e)
                if matches!(
                    e.kind(),
                    std::io::ErrorKind::WouldBlock | std::io::ErrorKind::TimedOut
                ) =>
            {
                if last_seen.elapsed() >= shared.cfg.idle_timeout {
                    break Ok(());
                }
            }
            Err(e) => break Err(e),
        }
    };

    alive.store(false, Ordering::Relaxed);
    let _ = pusher.join();
    shared.counters.agents.fetch_sub(1, Ordering::Relaxed);
    result
}

// ── client ──────────────────────────────────────────────────────────────────

/// A parsed `http://host:port[/prefix]` coordinator URL.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoordinatorUrl {
    /// `host:port`, as sent in the `Host:` header and dialled.
    pub host_port: String,
    /// Path prefix ahead of `/v1/…`, `""` when the hub is at the root.
    pub prefix: String,
}

impl CoordinatorUrl {
    /// Accepts `host:port`, `http://host:port`, and a path prefix for
    /// hubs published under a subpath by a reverse proxy.  `https://`
    /// is refused rather than silently downgraded: this client speaks
    /// no TLS, so an `https` URL must be reached through a local
    /// terminator (see `deploy/rho-coordinator/`).
    pub fn parse(url: &str) -> Result<Self, String> {
        let url = url.trim();
        if url.starts_with("https://") {
            return Err(format!(
                "{url}: this client speaks plain HTTP only; terminate TLS in front of it \
                 (stunnel/nginx) and pass the local http:// address"
            ));
        }
        let rest = url.strip_prefix("http://").unwrap_or(url);
        let rest = rest.trim_end_matches('/');
        let (host_port, prefix) = match rest.split_once('/') {
            Some((h, p)) => (h, format!("/{}", p.trim_matches('/'))),
            None => (rest, String::new()),
        };
        if host_port.is_empty() {
            return Err(format!("{url}: no host"));
        }
        let host_port = if host_port.contains(':') {
            host_port.to_string()
        } else {
            format!("{host_port}:80")
        };
        Ok(Self {
            host_port,
            prefix: if prefix == "/" { String::new() } else { prefix },
        })
    }

    fn route(&self, path: &str) -> String {
        format!("{}{}", self.prefix, path)
    }

    fn connect(&self, timeout: Duration) -> Result<TcpStream, String> {
        let addr = self
            .host_port
            .to_socket_addrs()
            .map_err(|e| format!("{}: {e}", self.host_port))?
            .next()
            .ok_or_else(|| format!("{}: no address", self.host_port))?;
        TcpStream::connect_timeout(&addr, timeout).map_err(|e| format!("{}: {e}", self.host_port))
    }
}

fn request(
    url: &CoordinatorUrl,
    token: Option<&str>,
    method: &str,
    path: &str,
    body: Option<&[u8]>,
) -> Result<(u16, String), String> {
    let mut stream = url.connect(Duration::from_secs(10))?;
    stream
        .set_read_timeout(Some(Duration::from_secs(60)))
        .map_err(|e| e.to_string())?;
    let mut head = format!(
        "{method} {} HTTP/1.1\r\nHost: {}\r\nConnection: close\r\n",
        url.route(path),
        url.host_port
    );
    if let Some(t) = token {
        head.push_str(&format!("Authorization: Bearer {t}\r\n"));
    }
    if let Some(b) = body {
        head.push_str(&format!(
            "Content-Type: application/json\r\nContent-Length: {}\r\n",
            b.len()
        ));
    }
    head.push_str("\r\n");
    stream
        .write_all(head.as_bytes())
        .map_err(|e| e.to_string())?;
    if let Some(b) = body {
        stream.write_all(b).map_err(|e| e.to_string())?;
    }
    stream.flush().map_err(|e| e.to_string())?;

    let mut reader = BufReader::new(stream);
    let mut status_line = String::new();
    reader
        .read_line(&mut status_line)
        .map_err(|e| e.to_string())?;
    let status: u16 = status_line
        .split_whitespace()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| format!("bad status line {status_line:?}"))?;
    let mut length: Option<usize> = None;
    loop {
        let mut line = String::new();
        let n = reader.read_line(&mut line).map_err(|e| e.to_string())?;
        if n == 0 || line == "\r\n" || line == "\n" {
            break;
        }
        if let Some((name, value)) = line.split_once(':') {
            if name.trim().eq_ignore_ascii_case("content-length") {
                length = value.trim().parse().ok();
            }
        }
    }
    let mut out = Vec::new();
    match length {
        Some(n) => {
            out.resize(n, 0);
            reader.read_exact(&mut out).map_err(|e| e.to_string())?;
        }
        None => {
            reader.read_to_end(&mut out).map_err(|e| e.to_string())?;
        }
    }
    String::from_utf8(out)
        .map(|s| (status, s))
        .map_err(|e| e.to_string())
}

/// Fetch the job document, so an agent can be started with a URL and
/// nothing else — the fleet image carries no job file.
pub fn fetch_job(url: &CoordinatorUrl, token: Option<&str>) -> Result<JobSpec, String> {
    let (status, body) = request(url, token, "GET", "/v1/job", None)?;
    if status != 200 {
        return Err(format!("GET /v1/job: HTTP {status}: {}", body.trim()));
    }
    JobSpec::from_json(&body)
}

/// Fetch `/v1/status`.
pub fn fetch_status(url: &CoordinatorUrl, token: Option<&str>) -> Result<HubStatus, String> {
    let (status, body) = request(url, token, "GET", "/v1/status", None)?;
    if status != 200 {
        return Err(format!("GET /v1/status: HTTP {status}: {}", body.trim()));
    }
    serde_json::from_str(&body).map_err(|e| e.to_string())
}

/// One `POST /v1/sync`: push what the hub lacks, merge what it returns.
/// This is the pollable fallback — [`ReverseChannel`] is the live path.
pub fn sync_once(
    url: &CoordinatorUrl,
    token: Option<&str>,
    ctx: &JobContext,
    state: &Mutex<SharedState>,
) -> Result<(usize, usize), String> {
    let (known, mine) = {
        let st = state.lock().unwrap();
        (st.version_vector(), st.delta_for(&BTreeMap::new()))
    };
    let mut received = 0;
    let mut rejected = 0;
    // Pull first: the hub's answer is what we lack.
    let pull = serde_json::to_vec(&Frame::Hello {
        job_id: ctx.job_id.clone(),
        peer: String::new(),
        known,
    })
    .map_err(|e| e.to_string())?;
    let (status, body) = request(url, token, "POST", "/v1/sync", Some(&pull))?;
    if status != 200 {
        return Err(format!("POST /v1/sync: HTTP {status}: {}", body.trim()));
    }
    match serde_json::from_str::<Frame>(&body).map_err(|e| e.to_string())? {
        Frame::Batch {
            checkins, known, ..
        } => {
            let mut st = state.lock().unwrap();
            for ci in &checkins {
                match st.apply(ctx, ci, now_secs(), true) {
                    Ok(o) => {
                        if o.new {
                            received += 1;
                        }
                        rejected += o.rejected_dps;
                    }
                    Err(_) => rejected += 1,
                }
            }
            // Push exactly what the hub is missing.
            let to_send = st.delta_for(&known);
            drop(st);
            if !to_send.is_empty() {
                let push = serde_json::to_vec(&Frame::Push { checkins: to_send })
                    .map_err(|e| e.to_string())?;
                let (status, body) = request(url, token, "POST", "/v1/sync", Some(&push))?;
                if status != 200 {
                    return Err(format!("POST /v1/sync: HTTP {status}: {}", body.trim()));
                }
            }
        }
        Frame::Error { message } => return Err(message),
        _ => {
            // A hub that answered something else is a version skew;
            // the unsent local delta is retried on the next sync.
            let _ = mine;
            return Err("unexpected reply to sync".into());
        }
    }
    Ok((received, rejected))
}

/// Agent-side counters for the reverse channel.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ChannelStats {
    /// Check-ins merged from the hub.
    pub received: u64,
    /// Check-ins sent up.
    pub sent: u64,
    /// Records the local state refused (should stay 0).
    pub rejected: u64,
    /// Times the channel was (re)established.
    pub connects: u64,
    /// Times it dropped and backed off.
    pub reconnects: u64,
    /// `true` while a socket is up.
    pub connected: bool,
    /// Last error seen, for the status line.
    pub last_error: Option<String>,
}

impl ChannelStats {
    fn note_error(&mut self, e: impl std::fmt::Display) {
        self.last_error = Some(e.to_string());
    }
}

/// The agent's end of the reverse channel: one outbound connection,
/// held open, reconnected with exponential backoff, feeding the local
/// [`SharedState`] and carrying this agent's check-ins up.
///
/// The agent needs no inbound reachability and no public address — it
/// dials the hub, and the hub answers on the same socket.
pub struct ReverseChannel {
    stats: Arc<Mutex<ChannelStats>>,
    outbox: Arc<Mutex<Vec<CheckIn>>>,
    stop: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
}

impl ReverseChannel {
    /// Start the background connector.  Returns immediately: the first
    /// connection (and every later one) happens on the worker thread,
    /// so a hub that is down or slow never blocks the walk.
    pub fn start(
        url: CoordinatorUrl,
        token: Option<String>,
        peer: String,
        ctx: Arc<JobContext>,
        state: Arc<Mutex<SharedState>>,
    ) -> Self {
        let stats = Arc::new(Mutex::new(ChannelStats::default()));
        let outbox = Arc::new(Mutex::new(Vec::new()));
        let stop = Arc::new(AtomicBool::new(false));
        let handle = {
            let (stats, outbox, stop) =
                (Arc::clone(&stats), Arc::clone(&outbox), Arc::clone(&stop));
            std::thread::spawn(move || {
                let mut backoff = Duration::from_secs(1);
                while !stop.load(Ordering::Relaxed) {
                    match run_agent_session(
                        &url,
                        token.as_deref(),
                        &peer,
                        &ctx,
                        &state,
                        &outbox,
                        &stats,
                        &stop,
                    ) {
                        Ok(()) => backoff = Duration::from_secs(1),
                        Err(e) => {
                            stats.lock().unwrap().note_error(e);
                            backoff = (backoff * 2).min(Duration::from_secs(30));
                        }
                    }
                    {
                        let mut s = stats.lock().unwrap();
                        s.connected = false;
                        s.reconnects += 1;
                    }
                    // Sleep in slices so a stop is prompt.
                    let deadline = Instant::now() + backoff;
                    while Instant::now() < deadline && !stop.load(Ordering::Relaxed) {
                        std::thread::sleep(Duration::from_millis(100));
                    }
                }
            })
        };
        Self {
            stats,
            outbox,
            stop,
            handle: Some(handle),
        }
    }

    /// Hand a check-in to the channel.  Non-blocking: it is queued and
    /// goes up on the next turn of the session loop, or on the next
    /// connection if the hub is unreachable right now.
    pub fn publish(&self, ci: &CheckIn) {
        self.outbox.lock().unwrap().push(ci.clone());
    }

    pub fn stats(&self) -> ChannelStats {
        self.stats.lock().unwrap().clone()
    }

    /// Block until the queue has drained up the channel or `timeout`
    /// passes; `true` if it drained.  Called before an agent exits so
    /// the last check-ins — the solution among them — are not lost
    /// with the process.
    pub fn flush(&self, timeout: Duration) -> bool {
        let deadline = Instant::now() + timeout;
        loop {
            if self.outbox.lock().unwrap().is_empty() {
                return true;
            }
            if Instant::now() >= deadline {
                return false;
            }
            std::thread::sleep(Duration::from_millis(50));
        }
    }

    pub fn stop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

impl Drop for ReverseChannel {
    fn drop(&mut self) {
        self.stop();
    }
}

/// One connection's lifetime, agent side.
#[allow(clippy::too_many_arguments)]
fn run_agent_session(
    url: &CoordinatorUrl,
    token: Option<&str>,
    peer: &str,
    ctx: &JobContext,
    state: &Mutex<SharedState>,
    outbox: &Mutex<Vec<CheckIn>>,
    stats: &Mutex<ChannelStats>,
    stop: &AtomicBool,
) -> Result<(), String> {
    let stream = url.connect(Duration::from_secs(10))?;
    stream
        .set_read_timeout(Some(Duration::from_millis(500)))
        .map_err(|e| e.to_string())?;
    let mut writer = stream.try_clone().map_err(|e| e.to_string())?;
    let mut head = format!(
        "GET {} HTTP/1.1\r\nHost: {}\r\nConnection: Upgrade\r\nUpgrade: {CHANNEL_PROTOCOL}\r\n",
        url.route("/v1/channel"),
        url.host_port
    );
    if let Some(t) = token {
        head.push_str(&format!("Authorization: Bearer {t}\r\n"));
    }
    head.push_str("\r\n");
    writer
        .write_all(head.as_bytes())
        .map_err(|e| e.to_string())?;
    writer.flush().map_err(|e| e.to_string())?;

    let mut reader = BufReader::new(stream);
    let mut status_line = String::new();
    reader
        .read_line(&mut status_line)
        .map_err(|e| e.to_string())?;
    if !status_line.contains(" 101 ") {
        // Drain the headers so the error message can carry the body.
        let mut rest = String::new();
        let _ = reader.read_to_string(&mut rest);
        return Err(format!(
            "channel refused: {}{}",
            status_line.trim(),
            if rest.trim().is_empty() {
                String::new()
            } else {
                format!(" — {}", rest.trim().lines().last().unwrap_or_default())
            }
        ));
    }
    loop {
        let mut line = String::new();
        let n = reader.read_line(&mut line).map_err(|e| e.to_string())?;
        if n == 0 || line == "\r\n" || line == "\n" {
            break;
        }
    }
    {
        let mut s = stats.lock().unwrap();
        s.connected = true;
        s.connects += 1;
        s.last_error = None;
    }

    let hello = Frame::Hello {
        job_id: ctx.job_id.clone(),
        peer: peer.to_string(),
        known: state.lock().unwrap().version_vector(),
    };
    write_frame(&mut writer, &hello).map_err(|e| e.to_string())?;

    let mut last_push = Instant::now();
    loop {
        if stop.load(Ordering::Relaxed) {
            return Ok(());
        }
        // Anything the lanes produced since the last turn goes up.
        let pending: Vec<CheckIn> = std::mem::take(&mut *outbox.lock().unwrap());
        if !pending.is_empty() {
            let n = pending.len() as u64;
            write_frame(&mut writer, &Frame::Push { checkins: pending })
                .map_err(|e| e.to_string())?;
            stats.lock().unwrap().sent += n;
            last_push = Instant::now();
        } else if last_push.elapsed() >= Duration::from_secs(20) {
            write_frame(&mut writer, &Frame::Ping).map_err(|e| e.to_string())?;
            last_push = Instant::now();
        }
        match read_frame(&mut reader) {
            Ok(None) => return Ok(()),
            Ok(Some(Frame::Batch { checkins, .. })) => {
                let mut st = state.lock().unwrap();
                let (mut got, mut bad) = (0u64, 0u64);
                for ci in &checkins {
                    match st.apply(ctx, ci, now_secs(), true) {
                        Ok(o) => {
                            if o.new {
                                got += 1;
                            }
                            bad += o.rejected_dps as u64;
                        }
                        Err(_) => bad += 1,
                    }
                }
                drop(st);
                let mut s = stats.lock().unwrap();
                s.received += got;
                s.rejected += bad;
            }
            Ok(Some(Frame::Error { message })) => return Err(message),
            Ok(Some(_)) => {}
            Err(e)
                if matches!(
                    e.kind(),
                    std::io::ErrorKind::WouldBlock | std::io::ErrorKind::TimedOut
                ) => {}
            Err(e) => return Err(e.to_string()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::job::{demo_curve, JobSpec};
    use super::super::worker::{run_lane, LaneOptions};
    use super::*;
    use num_bigint::BigUint;

    fn spec_for(secret: u32) -> JobSpec {
        let curve = demo_curve("demo-mid").unwrap();
        let q = curve
            .generator()
            .scalar_mul(&BigUint::from(secret), &curve.a_fe());
        let mut spec = JobSpec::new(&curve, &q, "hub", 9).unwrap();
        spec.dp_bits = 3;
        spec.unit_size = 16;
        spec
    }

    fn hub(spec: &JobSpec, token: Option<&str>) -> (Coordinator, Arc<JobContext>) {
        let ctx = Arc::new(spec.build().unwrap());
        let state = Arc::new(Mutex::new(SharedState::new(&ctx)));
        let cfg = CoordinatorConfig {
            token: token.map(str::to_string),
            push_interval: Duration::from_millis(20),
            ..CoordinatorConfig::default()
        };
        let c = Coordinator::start("127.0.0.1:0", Arc::clone(&ctx), spec, state, cfg, 120, None)
            .unwrap();
        (c, ctx)
    }

    #[test]
    fn url_parsing_covers_the_shapes_an_operator_types() {
        assert_eq!(
            CoordinatorUrl::parse("http://10.0.1.7:8080").unwrap(),
            CoordinatorUrl {
                host_port: "10.0.1.7:8080".into(),
                prefix: String::new()
            }
        );
        // Bare host:port, and a default port for a proxied hostname.
        assert_eq!(
            CoordinatorUrl::parse("hub.internal:8080")
                .unwrap()
                .host_port,
            "hub.internal:8080"
        );
        assert_eq!(
            CoordinatorUrl::parse("http://hub.example.com")
                .unwrap()
                .host_port,
            "hub.example.com:80"
        );
        // A hub published under a subpath keeps the prefix.
        let u = CoordinatorUrl::parse("http://hub:8080/rho/").unwrap();
        assert_eq!(u.prefix, "/rho");
        assert_eq!(u.route("/v1/channel"), "/rho/v1/channel");
        // https is refused loudly rather than downgraded silently.
        let e = CoordinatorUrl::parse("https://hub.example.com").unwrap_err();
        assert!(e.contains("plain HTTP only"), "{e}");
    }

    #[test]
    fn frame_round_trip() {
        let f = Frame::Hello {
            job_id: "abc".into(),
            peer: "agent.0".into(),
            known: BTreeMap::from([("agent.0".to_string(), 7u64)]),
        };
        let s = serde_json::to_string(&f).unwrap();
        assert!(s.contains("\"kind\":\"hello\""));
        assert_eq!(serde_json::from_str::<Frame>(&s).unwrap(), f);
    }

    #[test]
    fn an_agent_needs_only_the_url_to_learn_the_job() {
        let spec = spec_for(4_242);
        let (c, ctx) = hub(&spec, None);
        let url = CoordinatorUrl::parse(&c.url()).unwrap();
        let fetched = fetch_job(&url, None).unwrap();
        assert_eq!(fetched.build().unwrap().job_id, ctx.job_id);
        let st = fetch_status(&url, None).unwrap();
        assert_eq!(st.job_id, ctx.job_id);
        assert_eq!(st.progress.steps, 0);
    }

    #[test]
    fn a_bad_token_is_refused_and_health_still_answers() {
        let spec = spec_for(7);
        let (c, _) = hub(&spec, Some("s3cret"));
        let url = CoordinatorUrl::parse(&c.url()).unwrap();
        assert!(fetch_job(&url, None).unwrap_err().contains("401"));
        assert!(fetch_job(&url, Some("wrong")).unwrap_err().contains("401"));
        assert!(fetch_job(&url, Some("s3cret")).is_ok());
        // The load balancer's probe carries no credential.
        let (status, body) = request(&url, None, "GET", "/healthz", None).unwrap();
        assert_eq!(status, 200);
        assert!(body.contains("\"ok\":true"));
        assert!(c.stats().unauthorized >= 2);
    }

    #[test]
    fn a_channel_that_names_another_job_is_told_so() {
        let (c, _) = hub(&spec_for(11), None);
        let other = Arc::new(spec_for(12).build().unwrap());
        let state = Arc::new(Mutex::new(SharedState::new(&other)));
        let url = CoordinatorUrl::parse(&c.url()).unwrap();
        let outbox = Mutex::new(Vec::new());
        let stats = Mutex::new(ChannelStats::default());
        let stop = AtomicBool::new(false);
        let err = run_agent_session(&url, None, "x.0", &other, &state, &outbox, &stats, &stop)
            .unwrap_err();
        assert!(err.contains("different job"), "{err}");
    }

    /// The property the whole module exists for: two agents that never
    /// listen on anything, each dialling out to one hub, converge —
    /// including the solution, which reaches an agent that did not find
    /// it, over a socket that agent opened.
    #[test]
    fn two_dial_out_agents_converge_through_the_hub() {
        let secret = 88_888u32;
        let spec = spec_for(secret);
        let (c, ctx) = hub(&spec, Some("tok"));
        let url = CoordinatorUrl::parse(&c.url()).unwrap();

        let agents: Vec<Arc<Mutex<SharedState>>> = (0..2)
            .map(|_| Arc::new(Mutex::new(SharedState::new(&ctx))))
            .collect();
        let channels: Vec<ReverseChannel> = agents
            .iter()
            .enumerate()
            .map(|(i, st)| {
                ReverseChannel::start(
                    url.clone(),
                    Some("tok".into()),
                    format!("a{i}"),
                    Arc::clone(&ctx),
                    Arc::clone(st),
                )
            })
            .collect();

        let mut rounds = 0;
        loop {
            rounds += 1;
            for (i, st) in agents.iter().enumerate() {
                run_lane(
                    &ctx,
                    st,
                    &LaneOptions {
                        max_walkers: 8,
                        checkin_every: 4,
                        claim_window: 1,
                        ..LaneOptions::new(&format!("a{i}.0"))
                    },
                    &mut |ci| channels[i].publish(ci),
                    &|| false,
                );
            }
            std::thread::sleep(Duration::from_millis(60));
            if agents.iter().all(|s| s.lock().unwrap().solution.is_some()) {
                break;
            }
            assert!(rounds < 400, "no convergence through the hub");
        }
        for s in &agents {
            let st = s.lock().unwrap();
            assert_eq!(st.solution, Some(BigUint::from(secret)));
            assert_eq!(st.rejected_dps, 0);
        }
        let stats = c.stats();
        assert!(stats.accepted > 0, "hub merged nothing");
        assert!(
            stats.pushed > 0,
            "hub pushed nothing down the reverse channel"
        );
        assert!(channels
            .iter()
            .all(|ch| ch.stats().received > 0 || ch.stats().sent > 0));
    }

    /// The pollable fallback carries the same facts as the channel.
    #[test]
    fn one_shot_sync_moves_checkins_both_ways() {
        let spec = spec_for(31_337);
        let (c, ctx) = hub(&spec, None);
        let url = CoordinatorUrl::parse(&c.url()).unwrap();
        let a = Mutex::new(SharedState::new(&ctx));
        let b = Mutex::new(SharedState::new(&ctx));

        run_lane(
            &ctx,
            &a,
            &LaneOptions {
                max_walkers: 8,
                checkin_every: 4,
                ..LaneOptions::new("a.0")
            },
            &mut |_| {},
            &|| false,
        );
        sync_once(&url, None, &ctx, &a).unwrap();
        let (received, rejected) = sync_once(&url, None, &ctx, &b).unwrap();
        assert!(received > 0, "b learned nothing from the hub");
        assert_eq!(rejected, 0);
        assert_eq!(
            a.lock().unwrap().version_vector(),
            b.lock().unwrap().version_vector()
        );
    }

    /// Durability hook: what the hub accepts can be mirrored, so an
    /// instance replacement does not lose the log.
    #[test]
    fn accepted_checkins_reach_the_durability_hook() {
        let spec = spec_for(555);
        let ctx = Arc::new(spec.build().unwrap());
        let state = Arc::new(Mutex::new(SharedState::new(&ctx)));
        let seen = Arc::new(Mutex::new(Vec::<String>::new()));
        let sink = Arc::clone(&seen);
        let c = Coordinator::start(
            "127.0.0.1:0",
            Arc::clone(&ctx),
            &spec,
            state,
            CoordinatorConfig::default(),
            120,
            Some(Box::new(move |ci: &CheckIn| {
                sink.lock().unwrap().push(format!("{}#{}", ci.peer, ci.seq))
            })),
        )
        .unwrap();
        let url = CoordinatorUrl::parse(&c.url()).unwrap();
        let mine = Mutex::new(SharedState::new(&ctx));
        run_lane(
            &ctx,
            &mine,
            &LaneOptions {
                max_walkers: 4,
                checkin_every: 2,
                ..LaneOptions::new("a.0")
            },
            &mut |_| {},
            &|| false,
        );
        sync_once(&url, None, &ctx, &mine).unwrap();
        let seen = seen.lock().unwrap();
        assert!(!seen.is_empty(), "hook never ran");
        assert!(seen.iter().all(|s| s.starts_with("a.0#")));
    }
}

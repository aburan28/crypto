# Running the rho collaboration from one coordinator on EC2

**This package makes no performance claim.** Under `AGENTS.md` §3 it is
plumbing: no row on any scoreboard moves because of it. What it is
allowed to claim is the failure mode it removes, and that is exactly
one — *the fleet could not reach each other*.

The gossip transport in `src/cryptanalysis/pollard_collab/net.rs`
assumes every node can accept a connection. On a cloud fleet almost
none can: the workers are in private subnets, on spot instances, behind
NAT, in containers with no inbound rule. One host can be reached by
everybody — an instance with an elastic IP, or a load balancer. So the
topology that actually deploys is a hub, and the connections have to be
opened by the agents.

```
  agent (private subnet)  ────┐
  agent (spot fleet)      ────┼──►  https://rho.example.com   (EC2 + nginx)
  agent (laptop, NAT)     ────┘         GET /v1/channel
                                        Upgrade: rho-collab/1
                                   ◄─── the hub pushes back down
                                        the socket the agent opened
```

The last arrow is the point. Once an agent has dialled in, the hub uses
that same socket in the other direction — the **reverse channel** — to
push other agents' distinguished points and the solution to a machine
it could never have dialled. Nothing polls; an agent hears about a
collision found on another continent as soon as the hub has merged it.

## What the coordinator is not

It is a rendezvous, not an authority. It holds the same CRDT every
agent holds, verifies every distinguished point on arrival exactly as
an agent does (`apply(.., verify = true)`), and **hands out no work**:
agents still choose their own units under the leasing rules in
`state.rs`. Consequences worth being explicit about:

* **Losing the hub loses no work.** Agents keep walking their claimed
  units and keep their own DP tables; when the hub returns they
  reconverge, the way the design doc's §8 partition row describes. The
  hub is a single point of *reachability*, not of truth.
* **A restarted hub is not an empty hub**, as long as `--mailbox`
  points at a directory that survives the instance (an EBS volume, an
  EFS mount, anything that replicates a directory). Every accepted
  check-in is written there and the whole log is reloaded at start.
* **The trust model is unchanged.** The bearer token is access control,
  not integrity: it stops an unauthenticated stranger flooding the log,
  while every record in the log is still self-verifying, so a
  credentialled liar can still only waste their own time.

## The three commands

```bash
# 1. On any machine: write the job document.
crypto cryptanalysis rho-collab init --curve demo-40 --secret 1badc0de --out job.json

# 2. On the EC2 instance: the hub.
crypto cryptanalysis rho-collab coordinator \
    --job job.json --listen 127.0.0.1:8080 \
    --token-file /etc/rho/token --require-token \
    --mailbox /var/lib/rho/log

# 3. On every agent, anywhere: one URL is the whole configuration.
export RHO_COORDINATOR_URL=https://rho.example.com
export RHO_COORDINATOR_TOKEN=…
crypto cryptanalysis rho-collab work --node "$(hostname)" --threads "$(nproc)"

# And from anywhere with the token:
crypto cryptanalysis rho-collab status --coordinator "$RHO_COORDINATOR_URL"
```

Step 3 takes no `--job`, no `--listen` and no `--peer`: the agent
fetches the job document from the hub, so a fleet image is baked once
and aimed with an environment variable. `--coordinator` and `--token`
are flags too; the environment variables exist so a systemd unit and an
autoscaling group can carry them instead.

## Files here

| file | what it is |
|---|---|
| `rho-coordinator.service` | the hub unit: loopback bind, token file, durable mailbox, `Restart=always` |
| `rho-agent.service` | the agent unit: environment file, no inbound anything |
| `nginx.conf` | TLS in front, **with the upgrade headers the reverse channel needs** |
| `user-data.sh` | EC2 user-data: build, fetch job from S3 and token from Secrets Manager, start the unit |

## Wiring it up on AWS

**Security groups.** The hub instance needs inbound 443 from the agents
(or from the ALB's group) and nothing else. Agents need **no inbound
rule at all** — that is the whole reason for the reverse channel — and
outbound 443 to the hub.

**TLS.** The coordinator speaks plain HTTP on purpose: it is one
process, and a TLS stack inside it is a maintenance surface with no
upside when every deployment already has a terminator. Put nginx (see
`nginx.conf`) or an ALB in front and bind the hub to `127.0.0.1`. The
client refuses an `https://` URL rather than silently downgrading it,
so if agents are given one they must reach it through that terminator.

**Load balancers and proxies.** Two settings decide whether the reverse
channel works at all:

* the proxy must pass `Upgrade` and `Connection` through (`nginx.conf`
  does; an ALB does this for WebSocket-style upgrades natively);
* the idle timeout must exceed the hub's keepalive. The hub pings an
  idle channel every `idle_timeout/3` — 100 s at the default — so an
  ALB idle timeout of 300 s or more is safe. Below that, channels are
  closed under the agents, which costs a reconnect (with backoff) and
  nothing else, but fills the logs.

`/healthz` is deliberately outside the token so a target group can poll
it; it reveals nothing about the job.

**The token.** Generate it with `openssl rand -hex 32`, keep it in
Secrets Manager, and hand it to processes as a file or an environment
variable — never as a command-line argument, which `ps` shows to every
local user. `--require-token` makes the hub refuse to start without
one; set it whenever the bind address is reachable from outside the
host. A hub on a public address with no token warns loudly and keeps
going, because a private-subnet hub with a security group in front is a
legitimate configuration.

**Instance sizing.** The hub does two scalar multiplications per
distinguished point it accepts and holds the DP table in memory: at
`dp_bits` chosen by the default rule the table is about `n^{1/4}`
records of ~4 scalars, so a 2-vCPU instance carries a fleet whose
aggregate rate is a few thousand DPs a second. When the table stops
fitting, the answer is the sharding sketched in §10 of the design doc,
not a bigger instance.

**Spot agents.** An interrupted agent costs the unfinished part of its
current unit: its lease expires after `--lease-secs` and another agent
resumes the unit from the last cursor anybody reported. Nothing else is
needed to make the fleet spot-safe.

## Checking it end to end

```bash
crypto cryptanalysis rho-collab init --curve demo-40 --secret 1badc0de --out job.json
crypto cryptanalysis rho-collab coordinator --job job.json \
    --listen 127.0.0.1:18080 --token tok --mailbox ./log &
RHO_COORDINATOR_URL=http://127.0.0.1:18080 RHO_COORDINATOR_TOKEN=tok \
    crypto cryptanalysis rho-collab work --node alice --threads 2 &
RHO_COORDINATOR_URL=http://127.0.0.1:18080 RHO_COORDINATOR_TOKEN=tok \
    crypto cryptanalysis rho-collab work --node bob --threads 2
```

Both agents and the hub print the same `solution:` line with
`verified: true`, and neither agent ever listened on a port. On the
40-bit demo curve that takes about twenty seconds in a release build.

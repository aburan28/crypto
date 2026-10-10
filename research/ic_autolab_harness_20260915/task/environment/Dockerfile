FROM ubuntu:24.04@sha256:224a1869083a311ef3f13648a154ba79832fbef6364d31493642ca03082da254
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends build-essential pkg-config libssl-dev python3 valgrind git ca-certificates procps ripgrep && rm -rf /var/lib/apt/lists/*
COPY rust /opt/rust
COPY vendor /opt/vendor
COPY vendor-config.toml /opt/cargo/config.toml
COPY opencode /usr/local/bin/opencode
ENV PATH=/opt/rust/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV CARGO_HOME=/opt/cargo
COPY baseline /opt/baseline
COPY evaluator /opt/evaluator
RUN mkdir /opt/baseline-built && python3 -c "import sys; from pathlib import Path; sys.path.insert(0,'/opt/evaluator'); import tournament; tournament.snapshot_build(Path('/opt/baseline'),Path('/opt/baseline-built'))"
COPY harness /opt/harness
ENV CARGO_HOME=/home/researcher/.cargo
RUN useradd -m -u 1001 -s /bin/bash researcher && mkdir -p /app /home/researcher/.cargo && cp /opt/cargo/config.toml /home/researcher/.cargo/config.toml && cp -a /opt/baseline/. /app/ && cp /opt/harness/candidate.json /app/candidate.json && cp /opt/harness/AGENTS.md /app/AGENTS.md && chown -R researcher:researcher /app /home/researcher && chmod -R a+rX /opt && chmod a+rx /usr/local/bin/opencode
WORKDIR /app

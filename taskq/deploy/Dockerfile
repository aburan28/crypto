# taskq worker image. Toolchains for the experiments go in a derived image
# (FROM this one) per hardware class; this layer is only the queue client.
FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates build-essential \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/taskq
COPY pyproject.toml ./
COPY taskq ./taskq
RUN pip install --no-cache-dir '.[mcp]'
RUN useradd -m -u 10001 taskq && mkdir -p /var/cache/taskq && chown taskq /var/cache/taskq
USER taskq
ENTRYPOINT ["taskq"]
CMD ["worker", "--queue", "cpu"]

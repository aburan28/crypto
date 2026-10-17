FROM nvidia/cuda:13.3.1-runtime-ubuntu24.04

RUN rm -f /etc/apt/sources.list.d/cuda*.list \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
       ca-certificates curl libgomp1 python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/ecc2k130
COPY worker.py protocol.py seed_registry.py ./
COPY upload_gateway/bootstrap.py ./bootstrap.py

RUN chmod 0755 /opt/ecc2k130/bootstrap.py \
    && useradd --create-home --uid 10001 worker \
    && install -d -o worker -g worker /work/ecc2k130

USER worker
ENV ECC_ROOT=/work/ecc2k130 \
    PYTHONUNBUFFERED=1
ENTRYPOINT ["python3", "-u", "/opt/ecc2k130/bootstrap.py"]

# Build from repository root. One-shot process for an already authorized GPU pod.
FROM nvidia/cuda:12.8.1-devel-ubuntu22.04
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip libnuma1 procps util-linux && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir cupy-cuda12x==13.6.0 numpy==2.2.6
WORKDIR /study
COPY tools/curve_identity.py /study/tools/curve_identity.py
COPY tools/isolated_bench.py /study/tools/isolated_bench.py
COPY research/notes/ecc2k130/tau_adic_arithmetic_20260929 /study/research/notes/ecc2k130/tau_adic_arithmetic_20260929
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ENTRYPOINT ["python3", "research/notes/ecc2k130/tau_adic_arithmetic_20260929/gpu/launch_gpu.py"]
CMD ["--mode", "smoke", "--output", "/results/gpu-smoke.json"]

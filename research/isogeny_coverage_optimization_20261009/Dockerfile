FROM rust@sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends valgrind \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /work

# ── Stage 1: Build ────────────────────────────────────────────────────────────
FROM ubuntu:22.04 AS builder

ARG HTSLIB_VERSION=1.21
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -o Acquire::Check-Valid-Until=false update && apt-get install -y \
    curl \
    build-essential \
    pkg-config \
    clang \
    libclang-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libssl-dev \
    libdeflate-dev \
    && rm -rf /var/lib/apt/lists/*

# Build htslib from source at a pinned version
RUN curl -fsSL \
    https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    | tar -xj \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd .. \
    && rm -rf htslib-${HTSLIB_VERSION}

# Install Rust toolchain
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --default-toolchain stable --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

# Build the geac binary
WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src/ src/
RUN cargo build --release --jobs 4

# ── Stage 2: Runtime ──────────────────────────────────────────────────────────
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -o Acquire::Check-Valid-Until=false update && apt-get install -y \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    zlib1g \
    libssl3 \
    libdeflate0 \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Copy htslib shared libraries from builder
COPY --from=builder /usr/local/lib/libhts* /usr/local/lib/
RUN ldconfig

# Copy the geac binary
COPY --from=builder /build/target/release/geac /usr/local/bin/geac

# Install Python dependencies for the Explorer
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

# Copy the Explorer app
COPY app/ /app/

# Entrypoint script: "explorer" → streamlit, anything else → geac
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

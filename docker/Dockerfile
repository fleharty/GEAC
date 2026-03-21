# ── Stage 1: Build ────────────────────────────────────────────────────────────
# Pin to the BUILD platform so Rust compiles natively (no QEMU) even when the
# target is arm64.  Cross-compilation handles producing the right binary.
FROM --platform=$BUILDPLATFORM ubuntu:22.04 AS builder

ARG TARGETARCH
ARG HTSLIB_VERSION=1.21
ARG DEBIAN_FRONTEND=noninteractive

# Enable arm64 packages alongside native ones for cross-compilation deps
RUN dpkg --add-architecture arm64 || true

RUN apt-get -o Acquire::Check-Valid-Until=false update && apt-get install -y \
    curl \
    build-essential \
    pkg-config \
    clang \
    libclang-dev \
    # native (amd64) htslib build deps — also needed for bindgen header parsing
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libssl-dev \
    libdeflate-dev \
    # arm64 cross-compiler and arm64 htslib build deps
    gcc-aarch64-linux-gnu \
    g++-aarch64-linux-gnu \
    libbz2-dev:arm64 \
    liblzma-dev:arm64 \
    libcurl4-openssl-dev:arm64 \
    zlib1g-dev:arm64 \
    libssl-dev:arm64 \
    libdeflate-dev:arm64 \
    && rm -rf /var/lib/apt/lists/*

# Build htslib from source for the target architecture.
# Always installs to /usr/local so the runtime COPY path is unconditional.
RUN CC=$([ "$TARGETARCH" = "arm64" ] && echo "aarch64-linux-gnu-gcc" || echo "gcc") && \
    curl -fsSL \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    | tar -xj \
    && cd htslib-${HTSLIB_VERSION} \
    && CC=${CC} ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd .. \
    && rm -rf htslib-${HTSLIB_VERSION}

# Install Rust toolchain and add the arm64 target when cross-compiling
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --default-toolchain stable --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"
RUN [ "$TARGETARCH" = "arm64" ] && rustup target add aarch64-unknown-linux-gnu || true

# Build the geac binary
WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src/ src/
RUN if [ "$TARGETARCH" = "arm64" ]; then \
        CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER=aarch64-linux-gnu-gcc \
        PKG_CONFIG_ALLOW_CROSS=1 \
        PKG_CONFIG_PATH=/usr/local/lib/pkgconfig \
        HTSLIB_LIB_DIR=/usr/local/lib \
        HTSLIB_INCLUDE_DIR=/usr/local/include \
        BINDGEN_EXTRA_CLANG_ARGS="--target=aarch64-linux-gnu" \
        cargo build --release --target aarch64-unknown-linux-gnu --jobs 4 \
    && cp target/aarch64-unknown-linux-gnu/release/geac /geac; \
    else \
        cargo build --release --jobs 4 \
    && cp target/release/geac /geac; \
    fi

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

# Copy the geac binary (placed at /geac by the build step above)
COPY --from=builder /geac /usr/local/bin/geac

# Install Python dependencies for the Explorer
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

# Copy the Explorer app
COPY app/ /app/

# Entrypoint: "explorer" → streamlit, anything else → geac
RUN printf '#!/bin/bash\nif [ "$1" = "explorer" ]; then\n    shift\n    exec streamlit run /app/geac_explorer.py "$@"\nelse\n    exec geac "$@"\nfi\n' \
    > /usr/local/bin/entrypoint.sh && chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

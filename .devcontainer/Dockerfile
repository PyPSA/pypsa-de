FROM debian:bookworm

ARG PIXI_VERSION=v0.42.1
ENV XPAUTH_PATH=/mnt/storage/xpress/xpauth.xpr

RUN mkdir -p /mnt/storage

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates wget curl rsync git procps zsh build-essential python3-dev\
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Instal oh-my-zsh with git plugin and default theme
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.2.1/zsh-in-docker.sh)" -- \
    -t robbyrussell \
    -p git

RUN curl -L -o /usr/local/bin/pixi -fsSL --compressed "https://github.com/prefix-dev/pixi/releases/download/${PIXI_VERSION}/pixi-$(uname -m)-unknown-linux-musl" \
    && chmod +x /usr/local/bin/pixi \
    && pixi info
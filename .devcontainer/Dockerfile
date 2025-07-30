FROM debian:bookworm

ARG PIXI_VERSION=v0.49.0
ENV XPAUTH_PATH=/mnt/storage/xpress/xpauth.xpr

RUN mkdir -p /mnt/storage

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates wget curl rsync git delta procps zsh build-essential python3-dev\
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Instal oh-my-zsh with git plugin and default theme
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.2.1/zsh-in-docker.sh)" -- \
    -t robbyrussell \
    -p git

# Install Starship
RUN curl -sS https://starship.rs/install.sh | sh -s -- -y

RUN curl -L -o /usr/local/bin/pixi -fsSL --compressed "https://github.com/prefix-dev/pixi/releases/download/${PIXI_VERSION}/pixi-$(uname -m)-unknown-linux-musl" \
    && chmod +x /usr/local/bin/pixi \
    && pixi info

# build python environment: cannot find script
# RUN python3 /IdeaProjects/pypsa-at/scripts/pypsa-at/build_pixi_environment.py

# clone and install configuration files for git and other tools 
# RUN git clone https://gitlab.aggm.at/philip.worschischek/dotfiles.git ../dotfiles \
#     && ../dotfiles/setup.sh /.jbdevcontainer


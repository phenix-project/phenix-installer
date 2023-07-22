#!/bin/bash

set -xe

sudo mkdir -p /opt/empty_dir || true
for d in \
        /opt/ghc \
        /opt/hostedtoolcache \
        /usr/lib/jvm \
        /usr/local/.ghcup \
        /usr/local/android \
        /usr/local/powershell \
        /usr/share/dotnet \
        /usr/share/swift \
        ; do
  sudo rsync --stats -a --delete /opt/empty_dir/ $d || true
done
sudo apt-get purge -y -f firefox \
                          google-chrome-stable \
                          microsoft-edge-stable
sudo apt-get autoremove -y >& /dev/null
sudo apt-get autoclean -y >& /dev/null
sudo docker image prune --all --force
df -h

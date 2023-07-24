#!/bin/bash

# set -xe

sudo mkdir -p /opt/empty_dir || true
for d in \
        /opt/ghc \
        /opt/hostedtoolcache \
        /usr/lib/jvm \
        /usr/local/.ghcup \
        /usr/local/lib/android \
        /usr/local/share/powershell \
        /usr/share/dotnet \
        /usr/share/swift \
        ; do
  echo $d
  time sudo rsync --stats -a --delete /opt/empty_dir/ $d || true
  echo
done
time sudo apt-get purge -y -f firefox \
                          google-chrome-stable \
                          microsoft-edge-stable
time sudo apt-get autoremove -y >& /dev/null
time sudo apt-get autoclean -y >& /dev/null
echo
time sudo docker image prune --all --force
echo
df -h

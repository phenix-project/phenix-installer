#!/bin/bash

set -xe

# check disk space before
df -h
free -h

# clean image
mkdir -p /opt/empty_dir || true
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
  rsync --stats -a --delete /opt/empty_dir/ $d || true
done
apt-get purge -y -f firefox \
                    google-chrome-stable \
                    microsoft-edge-stable
apt-get autoremove -y >& /dev/null
apt-get autoclean -y >& /dev/null
docker image prune --all --force

rm -fr ${JAVA_HOME_8_X64}
rm -fr ${JAVA_HOME_11_X64}
rm -fr ${JAVA_HOME_17_X64}

rm -fr ${CHROMEWEBDRIVER}
rm -fr ${EDGEWEBDRIVER}
rm -fr ${GECKOWEBDRIVER}

# make swap disk
fallocate -l 16GiB /swapfile || true
chmod 600 /swapfile || true
mkswap /swapfile || true
swapon /swapfile

# check disk space after
df -h
free -h

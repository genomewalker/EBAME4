#!/usr/bin/env bash

# Install basic stuff
apt-get install cmake zlib1g-dev libbz2-dev

# Install MMSEQS2
cd /tmp || exit 'Folder not found'
git clone https://github.com/soedinglab/MMseqs2.git
cd MMseqs2 || exit 'Folder not found'
mkdir build
cd build || exit 'Folder not found'
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=/opt/MMseqs2 ..
make
sudo make install

echo 'export PATH=/opt/MMseqs2/bin/:$PATH' >> ~/.bashrc
echo 'source /opt/MMseqs2/util/bash-completion.sh' >> ~/.bashrc

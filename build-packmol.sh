#!/bin/bash

# Download file from GitHub and change name to packmol.tar.gz
curl -s https://api.github.com/repos/m3g/packmol/releases/latest | grep "tarball_url" | cut -d : -f 2,3 | tr -d \" | tr -d \, | wget -qi - -O packmol.tar.gz

# Unpack tarball and change directory name
tar -xzf packmol.tar.gz
mv m3g-packmol-* packmol

# Change directory and compile
cd packmol
make

# Make executable
chmod +x packmol

# Move executable to .local/bin
mv packmol ${HOME}/.local/bin

#!/bin/bash
set -e

# Enable nbgitpuller
jupyter serverextension enable --py nbgitpuller --sys-prefix


# Install LIGGGHTS
printf "Installing LIGGGHTS...\n"
current_dir=$(pwd)
parent_dir=${HOME}/PICI-LIGGGHTS
source_dir=$parent_dir/src
build_dir=$parent_dir/LIGGGHTS-build


# Create directories if needed
if [[ ! -d $parent_dir ]]; then
    git clone --depth 1 --branch 3.8.1 https://github.com/uob-positron-imaging-centre/PICI-LIGGGHTS.git $parent_dir
fi

if [[ ! -d $build_dir ]]; then
    mkdir $build_dir
fi


# Go into build dir and install LIGGGHTS
cd $build_dir
cmake -D CMAKE_BUILD_TYPE=Release $source_dir
make install -j 4
cd $current_dir

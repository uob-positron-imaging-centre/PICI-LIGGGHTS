#!/bin/bash
set -e

# Enable nbgitpuller
jupyter serverextension enable --py nbgitpuller --sys-prefix


# Install VTK
printf "Installing VTK...\n"
current_dir=$(pwd)
parent_dir=${HOME}
source_dir=$parent_dir/VTK-9.2.6
build_dir=$parent_dir/VTK-9.2.6/build


# Create directories if needed
if [[ ! -d $source_dir ]]; then
    git clone --depth 1 --branch v9.2.6 https://github.com/Kitware/VTK.git $source_dir
fi

if [[ ! -d $build_dir ]]; then
    mkdir $build_dir
fi


# Go into build dir and install VTK
cd $build_dir
cmake -D VTK_USE_MPI=ON -D VTK_GROUP_ENABLE_MPI=YES -D CMAKE_BUILD_TYPE=Release $source_dir
make install -j 10
cd $current_dir




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
make install -j 10
cd $current_dir

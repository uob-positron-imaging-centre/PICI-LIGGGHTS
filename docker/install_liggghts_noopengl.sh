#!/bin/bash
set -e

# Enable nbgitpuller
jupyter serverextension enable --py nbgitpuller --sys-prefix


# Install VTK
printf "Installing VTK...\n"
parentdir=$(pwd)
sourcedir=$parentdir/VTK-9.2.6
builddir=$parentdir/VTK-9.2.6-build


# Create directories if needed
if [[ ! -d $sourcedir ]]; then
    git clone --depth 1 --branch v9.2.6 https://github.com/Kitware/VTK.git $sourcedir
fi

if [[ ! -d $builddir ]]; then
    mkdir $builddir
fi


# Go into build dir and install VTK
cd $builddir
cmake \
    -D VTK_USE_MPI=ON                                           \
    -D VTK_GROUP_ENABLE_MPI=YES                                 \
    -D VTK_GROUP_ENABLE_Imaging=NO                              \
    -D VTK_GROUP_ENABLE_Qt=NO                                   \
    -D VTK_GROUP_ENABLE_StandAlone=NO                           \
    -D VTK_GROUP_ENABLE_Views=NO                                \
    -D VTK_GROUP_ENABLE_Web=NO                                  \
    -D CMAKE_BUILD_TYPE=Release                                 \
    -D VTK_MODULE_ENABLE_VTK_CommonCore=YES                     \
    -D VTK_MODULE_ENABLE_VTK_CommonExecutionModel=YES           \
    -D VTK_MODULE_ENABLE_VTK_CommonMath=YES                     \
    -D VTK_MODULE_ENABLE_VTK_CommonMisc=YES                     \
    -D VTK_MODULE_ENABLE_VTK_CommonSystem=YES                   \
    -D VTK_MODULE_ENABLE_VTK_CommonTransforms=YES               \
    -D VTK_MODULE_ENABLE_VTK_IOCore=YES                         \
    -D VTK_MODULE_ENABLE_VTK_IOLegacy=YES                       \
    -D VTK_MODULE_ENABLE_VTK_IOParallelXML=YES                  \
    -D VTK_MODULE_ENABLE_VTK_IOXML=YES                          \
    -D VTK_MODULE_ENABLE_VTK_IOXMLParser=YES                    \
    -D VTK_MODULE_ENABLE_VTK_ParallelCore=YES                   \
    -D VTK_MODULE_ENABLE_VTK_IONetCDF=YES                       \
    -D VTK_MODULE_ENABLE_VTK_IOGeometry=YES                     \
    -D VTK_MODULE_ENABLE_VTK_IOParallel=YES                     \
    -D VTK_MODULE_ENABLE_VTK_IOImage=YES                        \
    -D VTK_MODULE_ENABLE_VTK_FiltersCore=YES                    \
    -D VTK_MODULE_ENABLE_VTK_FiltersExtraction=YES              \
    -D VTK_MODULE_ENABLE_VTK_FiltersParallel=YES                \
    -D VTK_MODULE_ENABLE_VTK_FiltersGeneral=YES                 \
    -D VTK_MODULE_ENABLE_VTK_FiltersHybrid=YES                  \
    -D VTK_MODULE_ENABLE_VTK_ImagingCore=YES                    \
    -D VTK_MODULE_ENABLE_VTK_RenderingCore=YES                  \
    -D VTK_MODULE_ENABLE_VTK_FiltersVerdict=YES                 \
    -D VTK_MODULE_ENABLE_VTK_FiltersGeometry=YES                \
    -D VTK_MODULE_ENABLE_VTK_DomainsChemistry=YES               \
    -D VTK_MODULE_ENABLE_VTK_FiltersModeling=YES                \
    -D VTK_MODULE_ENABLE_VTK_FiltersSources=YES                 \
    -D VTK_MODULE_ENABLE_VTK_FiltersTexture=YES                 \
    -D VTK_MODULE_ENABLE_VTK_ImagingSources=YES                 \
    -D VTK_MODULE_ENABLE_VTK_FiltersStatistics=YES              \
    -D VTK_MODULE_ENABLE_VTK_ParallelDIY=YES                    \
    -D VTK_MODULE_ENABLE_VTK_CommonComputationalGeometry=YES    \
    $sourcedir

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

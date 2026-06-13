#!/bin/bash
###############################################################################
#
# Copyright 2008-2025, Mississippi State University
#
# This file is part of the Loci Framework.
#
# The Loci Framework is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The Loci Framework is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
#
###############################################################################

# Install with or without the Loci build information appended if set
# at configure time
if [ -z "${LOCI_INSTALL_DIR}" ]; then
    INSTALL_PATH=$INSTALL_DIR
else
    INSTALL_PATH=$INSTALL_DIR/$LOCI_INSTALL_DIR
fi

echo INSTALL_PATH = $INSTALL_PATH

# Copy files only if they differ or the destination does not exist. This is to
# avoid unnecessarily updating timestamps of files that have not changed,
# which can cause rebuilds of dependent files that have not changed.
soft_copy() {
    local src="$1"
    local dst="$2"

    mkdir -p "$(dirname "$dst")"
    if [ -e "$dst" ] && cmp -s "$src" "$dst"; then
        return 0
    fi

    cp -p "$src" "$dst"
}

# Copy a directory and its contents using soft_copy.
soft_copy_dir() {
    local src="$1"
    local dst_dir="$2"

    soft_copy "$src" "$dst_dir/$(basename "$src")"
}

# Copy files matching a glob pattern using soft_copy_dir.
glob_soft_copy() {
    local pattern="$1"
    local dst_dir="$2"
    local src

    while IFS= read -r src; do
        soft_copy_dir "$src" "$dst_dir"
    done < <(compgen -G "$pattern" || true)
}

echo Making Directories
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/lib
mkdir -p $INSTALL_PATH/bin

echo Installing Library Files
LIB_POSTFIX="so"

ARCH=${LOCI_ARCH-`uname -s`}
if [ $ARCH == "Darwin" ]; then
    LIB_POSTFIX="dylib"
fi
soft_copy_dir Tools/libTools.$LIB_POSTFIX $INSTALL_PATH/lib
soft_copy_dir System/libLoci.$LIB_POSTFIX $INSTALL_PATH/lib
soft_copy_dir FVMMod/fvm_m.so $INSTALL_PATH/lib
soft_copy_dir FVMMod/fvmFAD_m.so $INSTALL_PATH/lib
soft_copy_dir FVMMod/fvmVFAD_m.so $INSTALL_PATH/lib
soft_copy_dir FVMAdapt/fvmadapt_m.so $INSTALL_PATH/lib
soft_copy_dir FVMAdapt/libfvmadaptfunc.$LIB_POSTFIX $INSTALL_PATH/lib
soft_copy_dir FVMOverset/fvmoverset_m.so $INSTALL_PATH/lib
soft_copy_dir FVMOverset/fvmoversetFAD_m.so $INSTALL_PATH/lib
soft_copy_dir FVMOverset/fvmoversetVFAD_m.so $INSTALL_PATH/lib
soft_copy_dir sprng/libsprng.$LIB_POSTFIX $INSTALL_PATH/lib

if [ ! ${INSTALL_METIS} == 0 ]; then
    soft_copy_dir ParMetis-4.0/GKLib/libgk.$LIB_POSTFIX $INSTALL_PATH/lib
    soft_copy_dir ParMetis-4.0/METISLib/libmetis.$LIB_POSTFIX $INSTALL_PATH/lib
    soft_copy_dir ParMetis-4.0/ParMETISLib/libparmetis.$LIB_POSTFIX $INSTALL_PATH/lib
    mkdir -p $INSTALL_PATH/ParMetis-4.0
    mkdir -p $INSTALL_PATH/ParMetis-4.0/include
    mkdir -p $INSTALL_PATH/ParMetis-4.0/lib
    glob_soft_copy "ParMetis-4.0/*.h" $INSTALL_PATH/ParMetis-4.0
    glob_soft_copy "ParMetis-4.0/*.h" $INSTALL_PATH/ParMetis-4.0/include
    soft_copy_dir ParMetis-4.0/GKLib/libgk.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
    soft_copy_dir ParMetis-4.0/METISLib/libmetis.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
    soft_copy_dir ParMetis-4.0/ParMETISLib/libparmetis.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
fi

echo Installing Loci Tools
soft_copy_dir lpp/lpp $INSTALL_PATH/bin
soft_copy_dir FVMtools/cobalt2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/extract $INSTALL_PATH/bin
soft_copy_dir FVMtools/make_periodic $INSTALL_PATH/bin
soft_copy_dir FVMtools/plot3d2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/vog2surf $INSTALL_PATH/bin
soft_copy_dir FVMtools/ugrid2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/msh2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/foam2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/cfd++2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/fluent2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/ccm2vog $INSTALL_PATH/bin
soft_copy_dir FVMtools/vogmerge $INSTALL_PATH/bin
soft_copy_dir FVMtools/vogcheck $INSTALL_PATH/bin
soft_copy_dir FVMtools/vogcut $INSTALL_PATH/bin
soft_copy_dir FVMtools/extract_movie $INSTALL_PATH/bin
soft_copy_dir FVMtools/extruder $INSTALL_PATH/bin
soft_copy_dir FVMtools/refmesh $INSTALL_PATH/bin
soft_copy_dir FVMtools/marker $INSTALL_PATH/bin
soft_copy_dir FVMtools/refine $INSTALL_PATH/bin
soft_copy_dir FVMtools/cgns2ensight $INSTALL_PATH/bin
soft_copy_dir FVMtools/cgns2surf $INSTALL_PATH/bin
soft_copy_dir FVMtools/ugrid2cgns $INSTALL_PATH/bin
soft_copy_dir FVMtools/cgns2ugrid $INSTALL_PATH/bin
soft_copy_dir FVMtools/cgns2vog $INSTALL_PATH/bin

echo Copying config files: Loci.conf comp.conf sys.conf version.conf
soft_copy_dir Loci.conf $INSTALL_PATH
soft_copy_dir comp.conf $INSTALL_PATH
soft_copy_dir sys.conf $INSTALL_PATH
tmp_version_file=$(mktemp)
sed -e "s:^GIT_INFO.*:GIT_INFO = ${GIT_INFO}:" \
    -e "s:^GIT_BRANCH.*:GIT_BRANCH = ${GIT_BRANCH}:" \
		version.conf > $tmp_version_file
soft_copy $tmp_version_file $INSTALL_PATH/version.conf
rm -f $tmp_version_file

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
glob_soft_copy "include/*.h" $INSTALL_PATH/include
glob_soft_copy "include/*.lh" $INSTALL_PATH/include
soft_copy_dir include/Loci $INSTALL_PATH/include

for i in  Tools Config MPI_stubb FVMAdapt FVMOverset FVMMod; do
    mkdir -p $INSTALL_PATH/include/$i
    glob_soft_copy "include/$i/*.h" $INSTALL_PATH/include/$i
done
soft_copy_dir include/FVMOverset/overset $INSTALL_PATH/include/FVMOverset
glob_soft_copy "include/FVMOverset/*.lh" $INSTALL_PATH/include/FVMOverset
glob_soft_copy "include/FVMMod/*.lh" $INSTALL_PATH/include/FVMMod
glob_soft_copy "include/FVMAdapt/*.lh" $INSTALL_PATH/include/FVMAdapt

mkdir -p $INSTALL_PATH/docs
mkdir -p $INSTALL_PATH/docs/1D-Diffusion
mkdir -p $INSTALL_PATH/docs/heat
mkdir -p $INSTALL_PATH/docs/Datatypes
mkdir -p $INSTALL_PATH/docs/developer

if [ -e docs/tutorial/docs/tutorial.pdf ] ; then
    soft_copy_dir docs/tutorial/docs/tutorial.pdf $INSTALL_PATH/docs
fi
soft_copy_dir docs/tutorial/1D-Diffusion/Makefile $INSTALL_PATH/docs/1D-Diffusion
glob_soft_copy "docs/tutorial/1D-Diffusion/*.loci" $INSTALL_PATH/docs/1D-Diffusion
glob_soft_copy "docs/tutorial/1D-Diffusion/*.lh" $INSTALL_PATH/docs/1D-Diffusion
glob_soft_copy "docs/tutorial/1D-Diffusion/*.cc" $INSTALL_PATH/docs/1D-Diffusion
glob_soft_copy "docs/tutorial/1D-Diffusion/*.vars" $INSTALL_PATH/docs/1D-Diffusion

soft_copy_dir docs/tutorial/heat/Makefile $INSTALL_PATH/docs/heat
glob_soft_copy "docs/tutorial/heat/*.loci" $INSTALL_PATH/docs/heat
glob_soft_copy "docs/tutorial/heat/*.lh" $INSTALL_PATH/docs/heat
glob_soft_copy "docs/tutorial/heat/*.vog" $INSTALL_PATH/docs/heat
glob_soft_copy "docs/tutorial/heat/*.vars" $INSTALL_PATH/docs/heat

soft_copy_dir docs/tutorial/Datatypes/Makefile $INSTALL_PATH/docs/Datatypes
glob_soft_copy "docs/tutorial/Datatypes/*.cc" $INSTALL_PATH/docs/Datatypes

soft_copy_dir docs/developer/loci_style_guide.md $INSTALL_PATH/docs/developer

if [ -d docs/doxygen/html ] ; then
    mkdir -p $INSTALL_PATH/docs/doxygen
    rm -rf $INSTALL_PATH/docs/doxygen/html
    cp -a docs/doxygen/html $INSTALL_PATH/docs/doxygen/
fi

chmod -R a+rX $INSTALL_PATH

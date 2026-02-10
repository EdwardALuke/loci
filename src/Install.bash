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

copy_if_different() {
    local src="$1"
    local dst="$2"

    mkdir -p "$(dirname "$dst")"
    if [ -e "$dst" ] && cmp -s "$src" "$dst"; then
        return 0
    fi

    cp -p "$src" "$dst"
}

copy_into_dir_if_different() {
    local src="$1"
    local dst_dir="$2"

    copy_if_different "$src" "$dst_dir/$(basename "$src")"
}

copy_glob_into_dir_if_different() {
    local pattern="$1"
    local dst_dir="$2"
    local src

    while IFS= read -r src; do
        copy_into_dir_if_different "$src" "$dst_dir"
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
copy_into_dir_if_different Tools/libTools.$LIB_POSTFIX $INSTALL_PATH/lib
copy_into_dir_if_different System/libLoci.$LIB_POSTFIX $INSTALL_PATH/lib
copy_into_dir_if_different FVMMod/fvm_m.so $INSTALL_PATH/lib
copy_into_dir_if_different FVMMod/fvmFAD_m.so $INSTALL_PATH/lib
copy_into_dir_if_different FVMMod/fvmVFAD_m.so $INSTALL_PATH/lib
copy_into_dir_if_different FVMAdapt/fvmadapt_m.so $INSTALL_PATH/lib
copy_into_dir_if_different FVMAdapt/libfvmadaptfunc.$LIB_POSTFIX $INSTALL_PATH/lib
copy_into_dir_if_different FVMOverset/fvmoverset_m.so $INSTALL_PATH/lib
copy_into_dir_if_different FVMOverset/fvmoversetFAD_m.so $INSTALL_PATH/lib
copy_into_dir_if_different FVMOverset/fvmoversetVFAD_m.so $INSTALL_PATH/lib
copy_into_dir_if_different sprng/libsprng.$LIB_POSTFIX $INSTALL_PATH/lib

if [ ! ${INSTALL_METIS} == 0 ]; then
    copy_into_dir_if_different ParMetis-4.0/GKLib/libgk.$LIB_POSTFIX $INSTALL_PATH/lib
    copy_into_dir_if_different ParMetis-4.0/METISLib/libmetis.$LIB_POSTFIX $INSTALL_PATH/lib
    copy_into_dir_if_different ParMetis-4.0/ParMETISLib/libparmetis.$LIB_POSTFIX $INSTALL_PATH/lib
    mkdir -p $INSTALL_PATH/ParMetis-4.0
    mkdir -p $INSTALL_PATH/ParMetis-4.0/include
    mkdir -p $INSTALL_PATH/ParMetis-4.0/lib
    copy_glob_into_dir_if_different "ParMetis-4.0/*.h" $INSTALL_PATH/ParMetis-4.0
    copy_glob_into_dir_if_different "ParMetis-4.0/*.h" $INSTALL_PATH/ParMetis-4.0/include
    copy_into_dir_if_different ParMetis-4.0/GKLib/libgk.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
    copy_into_dir_if_different ParMetis-4.0/METISLib/libmetis.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
    copy_into_dir_if_different ParMetis-4.0/ParMETISLib/libparmetis.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
fi

echo Installing Loci Tools
copy_into_dir_if_different lpp/lpp $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/cobalt2vog $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/extract $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/make_periodic $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/plot3d2vog $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/vog2surf $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/ugrid2vog $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/msh2vog $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/cfd++2vog $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/fluent2vog $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/ccm2vog $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/vogmerge $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/vogcheck $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/vogcut $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/extract_movie $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/extruder $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/refmesh $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/marker $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/refine $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/cgns2ensight $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/cgns2surf $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/ugrid2cgns $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/cgns2ugrid $INSTALL_PATH/bin
copy_into_dir_if_different FVMtools/cgns2vog $INSTALL_PATH/bin

echo Copying config files: Loci.conf comp.conf sys.conf version.conf
copy_into_dir_if_different Loci.conf $INSTALL_PATH
copy_into_dir_if_different comp.conf $INSTALL_PATH
copy_into_dir_if_different sys.conf $INSTALL_PATH
tmp_version_file=$(mktemp)
sed -e "s:^GIT_INFO.*:GIT_INFO = ${GIT_INFO}:" \
    -e "s:^GIT_BRANCH.*:GIT_BRANCH = ${GIT_BRANCH}:" \
		version.conf > $tmp_version_file
copy_if_different $tmp_version_file $INSTALL_PATH/version.conf
rm -f $tmp_version_file

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
copy_glob_into_dir_if_different "include/*.h" $INSTALL_PATH/include
copy_glob_into_dir_if_different "include/*.lh" $INSTALL_PATH/include
copy_into_dir_if_different include/Loci $INSTALL_PATH/include

for i in  Tools Config MPI_stubb FVMAdapt FVMOverset FVMMod; do
    mkdir -p $INSTALL_PATH/include/$i
    copy_glob_into_dir_if_different "include/$i/*.h" $INSTALL_PATH/include/$i
done
copy_into_dir_if_different include/FVMOverset/overset $INSTALL_PATH/include/FVMOverset
copy_glob_into_dir_if_different "include/FVMOverset/*.lh" $INSTALL_PATH/include/FVMOverset
copy_glob_into_dir_if_different "include/FVMMod/*.lh" $INSTALL_PATH/include/FVMMod
copy_glob_into_dir_if_different "include/FVMAdapt/*.lh" $INSTALL_PATH/include/FVMAdapt

mkdir -p $INSTALL_PATH/docs
mkdir -p $INSTALL_PATH/docs/1D-Diffusion
mkdir -p $INSTALL_PATH/docs/heat
mkdir -p $INSTALL_PATH/docs/Datatypes

if [ -e Tutorial/docs/tutorial.pdf ] ; then
    copy_into_dir_if_different Tutorial/docs/tutorial.pdf $INSTALL_PATH/docs
fi
copy_into_dir_if_different Tutorial/1D-Diffusion/Makefile $INSTALL_PATH/docs/1D-Diffusion
copy_glob_into_dir_if_different "Tutorial/1D-Diffusion/*.loci" $INSTALL_PATH/docs/1D-Diffusion
copy_glob_into_dir_if_different "Tutorial/1D-Diffusion/*.lh" $INSTALL_PATH/docs/1D-Diffusion
copy_glob_into_dir_if_different "Tutorial/1D-Diffusion/*.cc" $INSTALL_PATH/docs/1D-Diffusion
copy_glob_into_dir_if_different "Tutorial/1D-Diffusion/*.vars" $INSTALL_PATH/docs/1D-Diffusion

copy_into_dir_if_different Tutorial/heat/Makefile $INSTALL_PATH/docs/heat
copy_glob_into_dir_if_different "Tutorial/heat/*.loci" $INSTALL_PATH/docs/heat
copy_glob_into_dir_if_different "Tutorial/heat/*.lh" $INSTALL_PATH/docs/heat
copy_glob_into_dir_if_different "Tutorial/heat/*.vog" $INSTALL_PATH/docs/heat
copy_glob_into_dir_if_different "Tutorial/heat/*.vars" $INSTALL_PATH/docs/heat

copy_into_dir_if_different Tutorial/Datatypes/Makefile $INSTALL_PATH/docs/Datatypes
copy_glob_into_dir_if_different "Tutorial/Datatypes/*.cc" $INSTALL_PATH/docs/Datatypes

chmod -R a+rX $INSTALL_PATH

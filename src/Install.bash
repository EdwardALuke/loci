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
cp Tools/libTools.$LIB_POSTFIX $INSTALL_PATH/lib
cp System/libLoci.$LIB_POSTFIX $INSTALL_PATH/lib
cp FVMMod/fvm_m.so $INSTALL_PATH/lib
cp FVMMod/fvmFAD_m.so $INSTALL_PATH/lib
cp FVMAdapt/fvmadapt_m.so $INSTALL_PATH/lib
cp FVMAdapt/libfvmadaptfunc.$LIB_POSTFIX $INSTALL_PATH/lib
cp FVMOverset/fvmoverset_m.so $INSTALL_PATH/lib
cp FVMOverset/fvmoversetFAD_m.so $INSTALL_PATH/lib
cp sprng/libsprng.$LIB_POSTFIX $INSTALL_PATH/lib

if [ ! ${INSTALL_METIS} == 0 ]; then
    cp ParMetis-4.0/GKLib/libgk.$LIB_POSTFIX $INSTALL_PATH/lib
    cp ParMetis-4.0/METISLib/libmetis.$LIB_POSTFIX $INSTALL_PATH/lib
    cp ParMetis-4.0/ParMETISLib/libparmetis.$LIB_POSTFIX $INSTALL_PATH/lib
    mkdir -p $INSTALL_PATH/ParMetis-4.0
    mkdir -p $INSTALL_PATH/ParMetis-4.0/include
    mkdir -p $INSTALL_PATH/ParMetis-4.0/lib
    cp ParMetis-4.0/*.h $INSTALL_PATH/ParMetis-4.0
    cp ParMetis-4.0/*.h $INSTALL_PATH/ParMetis-4.0/include
    cp ParMetis-4.0/GKLib/libgk.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
    cp ParMetis-4.0/METISLib/libmetis.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
    cp ParMetis-4.0/ParMETISLib/libparmetis.$LIB_POSTFIX $INSTALL_PATH/ParMetis-4.0/lib
fi

echo Installing Loci Tools
cp lpp/lpp $INSTALL_PATH/bin
cp FVMtools/cobalt2vog $INSTALL_PATH/bin
cp FVMtools/extract $INSTALL_PATH/bin
cp FVMtools/make_periodic $INSTALL_PATH/bin
cp FVMtools/plot3d2vog $INSTALL_PATH/bin
cp FVMtools/vog2surf $INSTALL_PATH/bin
cp FVMtools/ugrid2vog $INSTALL_PATH/bin
cp FVMtools/msh2vog $INSTALL_PATH/bin
cp FVMtools/cfd++2vog $INSTALL_PATH/bin
cp FVMtools/fluent2vog $INSTALL_PATH/bin
cp FVMtools/ccm2vog $INSTALL_PATH/bin
cp FVMtools/vogmerge $INSTALL_PATH/bin
cp FVMtools/vogcheck $INSTALL_PATH/bin
cp FVMtools/vogcut $INSTALL_PATH/bin
cp FVMtools/extract_movie $INSTALL_PATH/bin
cp FVMtools/extruder $INSTALL_PATH/bin
cp FVMtools/refmesh $INSTALL_PATH/bin
cp FVMtools/marker $INSTALL_PATH/bin
cp FVMtools/refine $INSTALL_PATH/bin
cp FVMtools/cgns2ensight $INSTALL_PATH/bin
cp FVMtools/cgns2surf $INSTALL_PATH/bin
cp FVMtools/ugrid2cgns $INSTALL_PATH/bin
cp FVMtools/cgns2ugrid $INSTALL_PATH/bin
cp FVMtools/cgns2vog $INSTALL_PATH/bin

echo Copying config files: Loci.conf comp.conf sys.conf version.conf
cp Loci.conf comp.conf sys.conf $INSTALL_PATH
sed -e "s:^GIT_INFO.*:GIT_INFO = ${GIT_INFO}:" \
    -e "s:^GIT_BRANCH.*:GIT_BRANCH = ${GIT_BRANCH}:" \
	version.conf > $INSTALL_PATH/version.conf

echo Installing \#include files
mkdir -p $INSTALL_PATH/include
cp include/*.h $INSTALL_PATH/include
cp include/*.lh $INSTALL_PATH/include
cp include/Loci $INSTALL_PATH/include

for i in  Tools Config MPI_stubb FVMAdapt FVMOverset; do
    mkdir -p $INSTALL_PATH/include/$i
    cp include/$i/*.h $INSTALL_PATH/include/$i
done
cp include/FVMOverset/*.lh $INSTALL_PATH/include/FVMOverset
cp include/FVMMod/*.lh $INSTALL_PATH/include/FVMMod

mkdir -p $INSTALL_PATH/docs
mkdir -p $INSTALL_PATH/docs/1D-Diffusion
mkdir -p $INSTALL_PATH/docs/heat
mkdir -p $INSTALL_PATH/docs/Datatypes

if [ -e Tutorial/docs/tutorial.pdf ] ; then
    cp Tutorial/docs/tutorial.pdf $INSTALL_PATH/docs
fi
cp Tutorial/1D-Diffusion/Makefile $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.loci $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.lh $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.cc $INSTALL_PATH/docs/1D-Diffusion
cp Tutorial/1D-Diffusion/*.vars $INSTALL_PATH/docs/1D-Diffusion

cp Tutorial/heat/Makefile $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.loci $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.lh $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.vog $INSTALL_PATH/docs/heat
cp Tutorial/heat/*.vars $INSTALL_PATH/docs/heat

cp Tutorial/Datatypes/Makefile $INSTALL_PATH/docs/Datatypes
cp Tutorial/Datatypes/*.cc $INSTALL_PATH/docs/Datatypes

chmod -R a+rX $INSTALL_PATH

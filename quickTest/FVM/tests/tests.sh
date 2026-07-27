#!/bin/bash
export LD_LIBRARY_PATH=$LOCI_BASE/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LOCI_BASE/lib:$DYLD_LIBRARY_PATH

BC_TESTS=$*

for ii in $BC_TESTS ; do
  i=`echo $ii | sed s/.vars//`
  cat $i.txt
  cat $i.txt > $i.test
  mkdir TEST_$i
  EXTRACT_ARGS=`cat $i.xtr.dat`
  TOL_ARGS=`cat $i.tol.dat`
  GRID_DAT=`cat $i.vog.dat`
  cd TEST_$i
  ln -s ../$i.vars .
  ln -s ../$GRID_DAT $i.vog
  if [ -f ../$i.sst.ic ]; then
    ln -s ../$i.sst.ic sst_hdf5.ic
  fi
  if [ -f ../$i.mdl ]; then
    ln -s ../$i.mdl $i.mdl
  fi
  $FLOWTEST $i > flow.out

  $EXTRACT -ascii $i 10ms $EXTRACT_ARGS > $i.dat
  cd ..
 if $NDIFF $i.dat TEST_$i/$i.dat $TOL_ARGS ; then
     echo Test Passed `pwd` $i
     echo Test `pwd` $i PASSED! >> $i.test
     rm -fr TEST_$i
 else 
     echo Test `pwd` $i FAILED !!!!!!!
     echo Test `pwd` $i FAILED! >> $i.test
  fi
done

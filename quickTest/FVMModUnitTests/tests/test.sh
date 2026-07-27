#!/bin/bash
export LD_LIBRARY_PATH=$LOCI_BASE/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LOCI_BASE/lib:$DYLD_LIBRARY_PATH
$MPIRUN ../FVMModUnitTest $LOCI_CMD_FLAGS $1 >& $2
tail -n 1 $2 > $2.tmp
if ! grep -q "PASSED" $2.tmp ; then echo TEST FVMModUnitTests/$1":" FAILED! >> $2; fi
rm $2.tmp


include $(LOCI_BASE)/Loci.conf
include $(LOCI_BASE)/version.conf
include $(TEST_BASE)/test.conf
LIB_PATH=$(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH=$(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)
FVMADAPT_TEST_BASE=$(TEST_BASE)/FVMAdaptTest
HEX_BLOCK=$(TEST_BASE)/FVMModUnitTests/tests/block.vog

default: testGridRef.dat testGridRef2.dat testGridRef3.dat testGridXml.dat hexBlockPar.dat
	cat dats/testGridRef.dat |sort -g > testGridRef.ref
	diff testGridRef.dat testGridRef.ref
	rm testGridRef.ref
	cat dats/testGridRef2.dat |sort -g > testGridRef2.ref
	diff testGridRef2.dat testGridRef2.ref
	rm testGridRef2.ref
	cat dats/testGridRef3.dat |sort -g > testGridRef3.ref
	diff testGridRef3.dat testGridRef3.ref
	rm testGridRef3.ref
	diff testGridXml.dat $(FVMADAPT_TEST_BASE)/dats/testGridXml.dat
	diff hexBlockPar.dat $(FVMADAPT_TEST_BASE)/dats/hexBlockPar.dat
	rm -fr output debug
	rm *.plan *.tag testGridRef*.vog *.quality
	rm testGridRef.dat testGridRef2.dat testGridRef3.dat
	rm testGridXml.vog testGridXml.log testGridXml.dat
	rm hexBlockPar.vog hexBlockPar.log hexBlockPar.dat

testGrid.tag: testGrid.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGrid
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(SERIALRUN) $(LOCI_BASE)/bin/extract -ascii testGrid 0 x | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/"| sed "s/.*e-.*/0/"> testGrid.tag

testGrid.plan: testGrid.vog testGrid.tag
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGrid.vog -tag testGrid.tag -o testGrid.plan

testGridRef.vog: testGrid.vog testGrid.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGrid.plan -o testGridRef.vog

output/testGridRef.topo: testGridRef.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGridRef

testGridRef.dat: output/testGridRef.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(SERIALRUN) $(LOCI_BASE)/bin/extract -ascii testGridRef 0 x y z | sort -g > testGridRef.dat

testGridRef.tag: output/testGridRef.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(SERIALRUN) $(LOCI_BASE)/bin/extract -ascii testGridRef 0 y | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/" |sed "s/.*e-.*/0/" > testGridRef.tag

testGridRef2.plan: testGrid.vog testGrid.plan testGridRef.tag
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGrid.plan -tag testGridRef.tag -o testGridRef2.plan

testGridRef2.vog: testGrid.vog testGridRef2.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGridRef2.plan -o testGridRef2.vog

output/testGridRef2.topo: testGridRef2.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGridRef2

testGridRef2.dat: output/testGridRef2.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(SERIALRUN) $(LOCI_BASE)/bin/extract -ascii testGridRef2 0 x y z | sort -g > testGridRef2.dat

testGridRef3.plan: testGridRef.vog testGridRef.tag
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGridRef.vog -tag testGridRef.tag -o testGridRef3.plan

testGridRef3.vog: testGridRef.vog testGridRef3.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGridRef.vog -r testGridRef3.plan -o testGridRef3.vog

output/testGridRef3.topo: testGridRef3.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGridRef3

testGridRef3.dat: output/testGridRef3.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(SERIALRUN) $(LOCI_BASE)/bin/extract -ascii testGridRef3 0 x y z | sort -g > testGridRef3.dat

testGridXml.plan: testGrid.vog $(FVMADAPT_TEST_BASE)/testGridRegion.xml
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGrid.vog \
	-xml $(FVMADAPT_TEST_BASE)/testGridRegion.xml -tol 0.08 -mode 2 -o testGridXml.plan

testGridXml.dat: testGridXml.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGrid.vog \
	-r testGridXml.plan -o testGridXml.vog > testGridXml.log 2>&1
	cat testGridXml.log
	grep -E '^num_(nodes|cells|faces) = ' testGridXml.log > testGridXml.dat

hexBlockPar.plan: $(HEX_BLOCK) $(FVMADAPT_TEST_BASE)/hexBlockSources.par
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g $(HEX_BLOCK) \
	-par $(FVMADAPT_TEST_BASE)/hexBlockSources.par -mode 2 -o hexBlockPar.plan

hexBlockPar.dat: hexBlockPar.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g $(HEX_BLOCK) \
	-r hexBlockPar.plan -o hexBlockPar.vog > hexBlockPar.log 2>&1
	cat hexBlockPar.log
	grep -E '^num_(nodes|cells|faces) = ' hexBlockPar.log > hexBlockPar.dat

clean:
	rm -fr output *.plan *.tag testGridRef*.vog *.quality testGridRef.dat testGridRef2.dat testGridRef3.dat *~ debug
	rm -f testGridXml.vog testGridXml.log testGridXml.dat
	rm -f hexBlockPar.vog hexBlockPar.log hexBlockPar.dat

distclean: clean


include $(TEST_BASE)/test.conf
LIB_PATH=$(LOCI_BASE)/lib:$(LD_LIBRARY_PATH)
DY_PATH=$(LOCI_BASE)/lib:$(DYLD_LIBRARY_PATH)

default: testGridRef.dat testGridRef2.dat testGridRef3.dat
	cat dats/testGridRef.dat |sort -g > testGridRef.ref
	diff testGridRef.dat testGridRef.ref
	rm testGridRef.ref
	cat dats/testGridRef2.dat |sort -g > testGridRef2.ref
	diff testGridRef2.dat testGridRef2.ref
	rm testGridRef2.ref
	cat dats/testGridRef3.dat |sort -g > testGridRef3.ref
	diff testGridRef3.dat testGridRef3.ref
	rm testGridRef3.ref
	rm -fr output debug
	rm *.plan *.tag testGridRef*.vog *.quality
	rm testGridRef.dat testGridRef2.dat testGridRef3.dat 

testGrid.tag: testGrid.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/vogcheck testGrid
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np 1 $(LOCI_BASE)/bin/extract -ascii testGrid 0 x | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/"| sed "s/.*e-.*/0/"> testGrid.tag

testGrid.plan: testGrid.vog testGrid.tag
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/marker -g testGrid.vog -tag testGrid.tag -o testGrid.plan

testGridRef.vog: testGrid.vog testGrid.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/refmesh -g testGrid.vog -r testGrid.plan -o testGridRef.vog

output/testGridRef.topo: testGridRef.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/vogcheck testGridRef

testGridRef.dat: output/testGridRef.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np 1 $(LOCI_BASE)/bin/extract -ascii testGridRef 0 x y z | sort -g > testGridRef.dat

testGridRef.tag: output/testGridRef.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np 1 $(LOCI_BASE)/bin/extract -ascii testGridRef 0 y | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/" |sed "s/.*e-.*/0/" > testGridRef.tag

testGridRef2.plan: testGrid.vog testGrid.plan testGridRef.tag
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/marker -g testGrid.vog -r testGrid.plan -tag testGridRef.tag -o testGridRef2.plan

testGridRef2.vog: testGrid.vog testGridRef2.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/refmesh -g testGrid.vog -r testGridRef2.plan -o testGridRef2.vog

output/testGridRef2.topo: testGridRef2.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/vogcheck testGridRef2

testGridRef2.dat: output/testGridRef2.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np 1 $(LOCI_BASE)/bin/extract -ascii testGridRef2 0 x y z | sort -g > testGridRef2.dat

testGridRef3.plan: testGridRef.vog testGridRef.tag
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/marker -g testGridRef.vog -tag testGridRef.tag -o testGridRef3.plan

testGridRef3.vog: testGridRef.vog testGridRef3.plan
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/refmesh -g testGridRef.vog -r testGridRef3.plan -o testGridRef3.vog

output/testGridRef3.topo: testGridRef3.vog
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np $(MPINUMPROC) $(LOCI_BASE)/bin/vogcheck testGridRef3

testGridRef3.dat: output/testGridRef3.topo
	LD_LIBRARY_PATH=$(LIB_PATH) DYLD_LIBRARY_PATH=$(DY_PATH) \
	$(MPIRUN) -np 1 $(LOCI_BASE)/bin/extract -ascii testGridRef3 0 x y z | sort -g > testGridRef3.dat

clean:
	rm -fr output *.plan *.tag testGridRef*.vog *.quality testGridRef.dat testGridRef2.dat testGridRef3.dat *~ debug

distclean: clean



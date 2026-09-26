CASES = restart
TOOLS = marker refmesh vogcheck extract
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
.NOTPARALLEL:
testGrid.vog:
	cp $(INPUT)/testGrid.vog $@

# The original regression: node-tag refinement, refinement from a saved plan,
# and refinement treating the generated mesh as a new base mesh.
restart: testGridRef.dat testGridRef2.dat testGridRef3.dat
	cat $(INPUT)/dats/testGridRef.dat |sort -g > testGridRef.ref
	diff testGridRef.dat testGridRef.ref
	cat $(INPUT)/dats/testGridRef2.dat |sort -g > testGridRef2.ref
	diff testGridRef2.dat testGridRef2.ref
	cat $(INPUT)/dats/testGridRef3.dat |sort -g > testGridRef3.ref
	diff testGridRef3.dat testGridRef3.ref

testGrid.tag: testGrid.vog
	$(PARALLEL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGrid
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii testGrid 0 x | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/"| sed "s/.*e-.*/0/"> testGrid.tag

testGrid.plan: testGrid.vog testGrid.tag
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGrid.vog -tag testGrid.tag -o testGrid.plan

testGridRef.vog: testGrid.vog testGrid.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGrid.plan -o testGridRef.vog

output/testGridRef.topo: testGridRef.vog
	$(PARALLEL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGridRef

testGridRef.dat: output/testGridRef.topo
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii testGridRef 0 x y z | sort -g > testGridRef.dat

testGridRef.tag: output/testGridRef.topo
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii testGridRef 0 y | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/" |sed "s/.*e-.*/0/" > testGridRef.tag

testGridRef2.plan: testGrid.vog testGrid.plan testGridRef.tag
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGrid.plan -tag testGridRef.tag -o testGridRef2.plan

testGridRef2.vog: testGrid.vog testGridRef2.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGridRef2.plan -o testGridRef2.vog

output/testGridRef2.topo: testGridRef2.vog
	$(PARALLEL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGridRef2

testGridRef2.dat: output/testGridRef2.topo
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii testGridRef2 0 x y z | sort -g > testGridRef2.dat

testGridRef3.plan: testGridRef.vog testGridRef.tag
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGridRef.vog -tag testGridRef.tag -o testGridRef3.plan

testGridRef3.vog: testGridRef.vog testGridRef3.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGridRef.vog -r testGridRef3.plan -o testGridRef3.vog

output/testGridRef3.topo: testGridRef3.vog
	$(PARALLEL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) testGridRef3

testGridRef3.dat: output/testGridRef3.topo
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii testGridRef3 0 x y z | sort -g > testGridRef3.dat

endif

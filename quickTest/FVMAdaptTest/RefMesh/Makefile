CASES = restart
TOOLS = marker refmesh vogcheck extract
NEEDS_H5DUMP = 1
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
.NOTPARALLEL:
.SECONDARY:
testGrid.vog:
	cp $(INPUT)/testGrid.vog $@

# Node tags refine the original mesh, its saved history, and its emitted VOG.
restart: testGrid.dat testGrid.summary testGridRef.dat testGridRef.summary \
         testGridRef2.dat testGridRef2.summary testGridRef3.dat testGridRef3.summary
	sh $(INPUT)/check_mesh.sh testGrid testGridRef
	sh $(INPUT)/check_mesh.sh testGridRef testGridRef2
	sh $(INPUT)/check_mesh.sh testGridRef testGridRef3

output/%.topo: %.vog
	$(PARALLEL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) $* > $*.check 2>&1

%.dat: output/%.topo
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii $* 0 x y z | sort -g > $@

%.summary: output/%.topo
	$(MESH_SUMMARY) $*.vog $*.check > $@

testGrid.tag: output/testGrid.topo
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii testGrid 0 x | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/" | sed "s/.*e-.*/0/" > $@

testGrid.plan: testGrid.vog testGrid.tag
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGrid.vog -tag testGrid.tag -mode 2 -o $@

testGridRef.vog: testGrid.vog testGrid.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGrid.plan -o $@

testGridRef.tag: output/testGridRef.topo
	$(SERIAL) $(LOCI_BASE)/bin/extract -ascii testGridRef 0 y | sed "s/0[.][01234].*/0/" | sed "s/0[.][56789].*/1/" | sed "s/.*e-.*/0/" > $@

testGridRef2.plan: testGrid.vog testGrid.plan testGridRef.tag
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGrid.plan -tag testGridRef.tag -mode 2 -o $@

testGridRef2.vog: testGrid.vog testGridRef2.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGrid.vog -r testGridRef2.plan -o $@

testGridRef3.plan: testGridRef.vog testGridRef.tag
	$(PARALLEL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g testGridRef.vog -tag testGridRef.tag -mode 2 -o $@

testGridRef3.vog: testGridRef.vog testGridRef3.plan
	$(PARALLEL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g testGridRef.vog -r testGridRef3.plan -o $@
endif

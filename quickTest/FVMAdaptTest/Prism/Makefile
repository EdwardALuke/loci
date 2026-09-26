CASES = round_trip empty_ranks
TOOLS = ugrid2vog marker refmesh vogcheck
NEEDS_H5DUMP = 1
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
prism.vog:
	cp $(INPUT)/prism.ugrid .
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m prism

# Refine one prism, then restore the base mesh by derefining its saved plan.
round_trip: prism.vog
	awk 'BEGIN { for (i=0; i<6; ++i) print -1 }' > derefine.tag
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g prism.vog \
	  -ctag $(INPUT)/refine.tag -mode 2 -o refined.plan
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g prism.vog -r refined.plan -o refined.vog
	$(SERIAL) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g prism.vog \
	  -r refined.plan -ctag derefine.tag -mode 2 -o restored.plan
	$(SERIAL) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g prism.vog -r restored.plan -o restored.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) prism > prism.check 2>&1
	$(MESH_SUMMARY) prism.vog prism.check > prism.actual
	diff -u $(INPUT)/dats/prism.dat prism.actual
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) refined > refined.check 2>&1
	$(MESH_SUMMARY) refined.vog refined.check > refined.actual
	diff -u $(INPUT)/dats/prismRef.dat refined.actual
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) restored > restored.check 2>&1
	$(MESH_SUMMARY) restored.vog restored.check > restored.actual
	diff -u $(INPUT)/dats/prismRoundTrip.dat restored.actual

# More ranks than cells must still produce the same refined mesh as serial.
empty_ranks: prism.vog
	$(THREE_RANKS) $(LOCI_BASE)/bin/marker $(LOCI_CMD_FLAGS) -g prism.vog \
	  -ctag $(INPUT)/refine.tag -mode 2 -o refined.plan
	$(THREE_RANKS) $(LOCI_BASE)/bin/refmesh $(LOCI_CMD_FLAGS) -g prism.vog -r refined.plan -o refined.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) refined > refined.check 2>&1
	$(MESH_SUMMARY) refined.vog refined.check > refined.actual
	diff -u $(INPUT)/dats/prismRef.dat refined.actual
endif

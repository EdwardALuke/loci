CASES = prism hex
TOOLS = ugrid2vog refine vogcheck
NEEDS_H5DUMP = 1
include $(dir $(lastword $(MAKEFILE_LIST)))../common.mk

ifdef CASE_WORK
REF_prism = mixedRefPrism.dat
REF_hex = mixedRefHex.dat

mixed.vog:
	cp $(INPUT)/mixed.ugrid .
	$(SERIAL) $(LOCI_BASE)/bin/ugrid2vog $(LOCI_CMD_FLAGS) -o -m mixed

# Refine either cell independently across the shared hex/prism face.
$(CASES): mixed.vog
	$(SERIAL) $(LOCI_BASE)/bin/refine $(LOCI_CMD_FLAGS) -g mixed.vog \
	  -ctag $(INPUT)/refine_$@.tag -mode 2 -o refined.vog
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) mixed > mixed.check 2>&1
	$(MESH_SUMMARY) mixed.vog mixed.check > mixed.actual
	diff -u $(INPUT)/dats/mixed.dat mixed.actual
	$(SERIAL) $(LOCI_BASE)/bin/vogcheck $(LOCI_CMD_FLAGS) refined > refined.check 2>&1
	$(MESH_SUMMARY) refined.vog refined.check > refined.actual
	diff -u $(INPUT)/dats/$(REF_$@) refined.actual
endif

TEST_BASE ?= $(abspath ..)
TESTS = RefMesh PlanInputs UniformHex Prism Mixed ExtrudedModes Core FaceOutput CellParents
.DEFAULT_GOAL := TestResults
.PHONY: TestResults $(TESTS) clean distclean

# Collect a summary even when an independent group fails.
TestResults:
	@rm -f $@
	@status=0; for group in $(TESTS); do \
	  rm -f $$group/TestResults; \
	  $(MAKE) --no-print-directory -C $$group LOCI_BASE="$(LOCI_BASE)" TEST_BASE="$(TEST_BASE)" || status=1; \
	  if test -f $$group/TestResults; then cat $$group/TestResults >> $@; \
	  else echo "FVMAdaptTest/$$group: FAILED (setup)" >> $@; fi; \
	done; printf '\nFVMAdapt summary\n'; cat $@; exit $$status

$(TESTS):
	@rm -f TestResults
	@$(MAKE) --no-print-directory -C $@ LOCI_BASE="$(LOCI_BASE)" TEST_BASE="$(TEST_BASE)"

# No compiler, MPI launcher, or Loci configuration is needed to clean.
clean distclean:
	rm -f TestResults
	@for group in $(TESTS); do $(MAKE) --no-print-directory -C $$group clean || exit $$?; done

TEST_BASE ?= $(abspath $(dir $(firstword $(MAKEFILE_LIST)))/..)
TESTS = RefMesh PlanInputs UniformHex Prism Mixed ExtrudedModes Core FaceOutput CellParents
GROUP_RESULTS := $(foreach t,$(TESTS),$t/TestResults)
CLEAN_TESTS := $(addprefix clean-,$(TESTS))
.DEFAULT_GOAL := TestResults
.PHONY: TestResults $(TESTS) $(CLEAN_TESTS) clean distclean FORCE

# Collect a summary even when an independent group fails.
TestResults: $(GROUP_RESULTS)
	@rm -f $@
	@cat $^ > $@
	@! grep -q FAILED $@

$(TESTS): %: %/TestResults
	@! grep -q FAILED $<

$(GROUP_RESULTS): %/TestResults: FORCE
	@mkdir -p $*
	@rm -f $@
	@$(MAKE) --no-print-directory -C $* -f "$(TEST_BASE)/FVMAdaptTest/$*/Makefile" LOCI_BASE="$(LOCI_BASE)" TEST_BASE="$(TEST_BASE)" \
	  || { test -s $@ || echo "FVMAdaptTest/$*: FAILED (make error)" > $@; }

FORCE:

# No compiler, MPI launcher, or Loci configuration is needed to clean.
clean distclean: $(CLEAN_TESTS)
	rm -f TestResults

$(CLEAN_TESTS): clean-%:
	@if test -d $*; then $(MAKE) --no-print-directory -C $* -f "$(TEST_BASE)/FVMAdaptTest/$*/Makefile" clean; fi

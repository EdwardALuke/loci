#TEST_BASE ?= $(abspath ..)
TESTS = RefMesh PlanInputs UniformHex Prism Mixed ExtrudedModes Core FaceOutput CellParents
GROUP_RESULTS := $(foreach t,$(TESTS),$t/TestResults)
CLEAN_TESTS := $(addprefix clean-,$(TESTS))
.DEFAULT_GOAL := TestResults
.PHONY: TestResults $(TESTS) $(CLEAN_TESTS) clean distclean FORCE

# Collect a summary even when an independent group fails.
TestResults: $(GROUP_RESULTS)
	@rm -f $@
	@cat $^ > $@
	@grep PASS $@
	@grep FAIL $@ || true
	@! grep -q FAIL $@

$(TESTS): %: %/TestResults
	@:

$(GROUP_RESULTS): %/TestResults: FORCE
	@$(MAKE) --no-print-directory -C $* LOCI_BASE="$(LOCI_BASE)" TEST_BASE="$(TEST_BASE)"

FORCE:

# No compiler, MPI launcher, or Loci configuration is needed to clean.
clean distclean: $(CLEAN_TESTS)
	rm -f TestResults

$(CLEAN_TESTS): clean-%:
	@$(MAKE) --no-print-directory -C $* clean

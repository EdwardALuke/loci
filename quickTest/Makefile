TEST_BASE = $(shell pwd)

.PHONY: FRC FVMModUnitTests FVMAdaptTests FVMAdaptTest FVMAdaptIllustrations Containers Tools SystemTests test default

default: FVM FVMModUnitTests FVMAdaptTests Containers Tools SystemTests
	cat Containers/TestResults Tools/TestResults FVMModUnitTests/TestResults FVM/TestResults FVMAdaptTests/TestResults SystemTests/TestResults > TestResults
	rm Containers/TestResults Tools/TestResults FVMModUnitTests/TestResults FVM/TestResults FVMAdaptTests/TestResults SystemTests/TestResults
	@grep PASS TestResults
	@grep FAIL TestResults || true
	@! grep -q FAIL TestResults

test: default

FVMModUnitTests: FRC
	$(MAKE) -C FVMModUnitTests tests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

FVM: FRC
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

FVMAdaptTests: FRC
	$(MAKE) -C FVMAdaptTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

FVMAdaptTest: FVMAdaptTests

FVMAdaptIllustrations: FRC
	$(MAKE) -C FVMAdaptTests Illustrations LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

Containers: FRC
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

Tools: FRC
	$(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

SystemTests: FRC
	$(MAKE) -C SystemTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

clean: FRC
	rm -f TestResults
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C FVMAdaptTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C SystemTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean

distclean:
	rm -f TestResults
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)  distclean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C FVMAdaptTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C SystemTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean

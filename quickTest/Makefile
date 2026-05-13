TEST_BASE = $(shell pwd)
default: FVM FVMModUnitTests FVMAdaptTest Containers Tools
	cat Containers/TestResults Tools/TestResults FVMModUnitTests/TestResults FVM/TestResults > TestResults
	rm Containers/TestResults Tools/TestResults FVMModUnitTests/TestResults FVM/TestResults
	@grep PASS TestResults
	@grep FAIL TestResults || /bin/true
	@! grep -q FAIL TestResults


.PHONY: FRC

FVMModUnitTests: FRC
	$(MAKE) -C FVMModUnitTests tests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

FVM: FRC
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

FVMAdaptTest: FRC
	$(MAKE) -C FVMAdaptTest LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

Containers: FRC
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

Tools: FRC
	$(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

clean: FRC
	rm -f TestResults
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C FVMAdaptTest LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean

distclean:
	rm -f TestResults
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)  distclean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C FVMAdaptTest LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean

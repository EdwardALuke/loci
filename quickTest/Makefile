TEST_BASE = $(shell pwd)
QUICKTEST_SRC = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
default: FVM FVMModUnitTests FVMAdaptTest Containers Tools

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
	@if [ ! -d Tools ]; then \
	  mkdir -p Tools ; \
	  ln -sf $(QUICKTEST_SRC)/Tools/Makefile Tools/Makefile ; \
	fi
	$(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)

clean: FRC
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean
	@if [ -d Tools ]; then \
	  $(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean ; \
	fi
	$(MAKE) -C FVMAdaptTest LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) clean

distclean:
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE)  distclean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean
	@if [ -d Tools ]; then \
	  $(MAKE) -C Tools LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean ; \
	fi
	$(MAKE) -C FVMAdaptTest LOCI_BASE=$(LOCI_BASE) TEST_BASE=$(TEST_BASE) distclean

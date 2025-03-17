default: FVM FVMModUnitTests Containers

.PHONY: FRC

FVMModUnitTests: FRC
	$(MAKE) -C FVMModUnitTests tests LOCI_BASE=$(LOCI_BASE)

FVM: FRC
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE)

Containers: FRC
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE)

clean: FRC
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) clean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) clean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) clean

distclean:
	$(MAKE) -C FVM LOCI_BASE=$(LOCI_BASE) distclean
	$(MAKE) -C FVMModUnitTests LOCI_BASE=$(LOCI_BASE) distclean
	$(MAKE) -C Containers LOCI_BASE=$(LOCI_BASE) distclean


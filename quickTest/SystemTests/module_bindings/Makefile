include $(LOCI_BASE)/Loci.conf

INCLUDES = -I../../contrib/doctest
MODULES = binding_consumer_m.so binding_provider_m.so binding_writer_m.so

.PHONY: all check check-mpi clean
all: test_bindings $(MODULES)

test_bindings: test_bindings.o
	$(LD) -o $@ $^ $(LIBS) $(LDFLAGS)

binding_%_m.so: %_lo.o
	$(SHARED_LD) $(SHARED_LD_FLAGS) $@ $(LIB_FLAGS) $^ $(BASE_LIBS)

# Keep debug/schedule output per process count. A failed assertion fails make.
check: all
	mkdir -p work/serial
	cd work/serial && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  ../../test_bindings --scheduleoutput

check-mpi: all
	mkdir -p work/mpi2
	cd work/mpi2 && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  $(MPI_RUN) -np 2 ../../test_bindings --scheduleoutput

clean:
	rm -f *.o *.d *.lcc~ $(MODULES) test_bindings
	rm -rf work

-include test_bindings.d consumer.d provider.d writer.d

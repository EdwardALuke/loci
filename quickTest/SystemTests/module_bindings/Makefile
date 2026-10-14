include $(LOCI_BASE)/Loci.conf

INCLUDES = -I../../contrib/doctest
MODULES = binding_consumer_m.so binding_provider_m.so binding_writer_m.so \
          binding_dummy_fluid_m.so binding_dummy_solid_m.so \
          binding_canonical_writer_m.so binding_dummy_field_m.so

.PHONY: all check check-mpi check-loading clean
all: test_bindings test_composition $(MODULES)

test_bindings: test_bindings.o
	$(LD) -o $@ $^ $(LIBS) $(LDFLAGS)

test_composition: test_composition.o
	$(LD) -o $@ $^ $(LIBS) $(LDFLAGS)

binding_%_m.so: %_lo.o
	$(SHARED_LD) $(SHARED_LD_FLAGS) $@ $(LIB_FLAGS) $^ $(BASE_LIBS)

# Keep debug/schedule output per process count. A failed assertion fails make.
check: all check-loading
	mkdir -p work/serial
	cd work/serial && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  ../../test_bindings --scheduleoutput
	cd work/serial && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  ../../test_composition --scheduleoutput
	mkdir -p work/loop-serial
	cd work/loop-serial && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  ../../test_composition --scheduleoutput --test-case="host exchange*"
	grep -Fq 'Iteration Loop{k}' work/loop-serial/debug/schedule
	grep -Fq 'Copy(fluid@wallTemperature){k}<-exchange@wallTemperature{k}' \
	  work/loop-serial/debug/schedule

# Each legacy-first mode needs a fresh process and an initially empty cache.
check-loading: all
	@set -e; for mode in rule-only full; do \
	  mkdir -p work/load-$$mode; \
	  (cd work/load-$$mode && \
	    LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	    ../../test_composition --scheduleoutput \
	      --test-case="loading after*" --subcase="$$mode*"); \
	done

check-mpi: all
	mkdir -p work/mpi2
	cd work/mpi2 && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  $(MPI_RUN) -np 2 ../../test_bindings --scheduleoutput
	cd work/mpi2 && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  $(MPI_RUN) -np 2 ../../test_composition --scheduleoutput
	mkdir -p work/loop-mpi2 work/map-mpi2
	cd work/loop-mpi2 && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  $(MPI_RUN) -np 2 ../../test_composition --scheduleoutput \
	    --test-case="host exchange*"
	cd work/map-mpi2 && \
	  LD_LIBRARY_PATH=$(LOCI_BASE)/lib:../..:$(LD_LIBRARY_PATH) \
	  $(MPI_RUN) -np 2 ../../test_composition --scheduleoutput \
	    --test-case="bound maps*"
	grep -Fq 'host@sample  to 1' work/map-mpi2/debug/schedule-0
	grep -Fq 'host@sample  from 0' work/map-mpi2/debug/schedule-1

clean:
	rm -f *.o *.d *.lcc~ $(MODULES) test_bindings test_composition
	rm -rf work

-include test_bindings.d consumer.d provider.d writer.d test_composition.d \
         dummy_fluid.d dummy_solid.d canonical_writer.d dummy_field.d

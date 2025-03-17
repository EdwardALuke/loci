include $(LOCI_BASE)/Loci.conf

TARGET = libsprng.$(LIB_SUFFIX)

INCLUDES=-I./include -DSPRNG_MPI -DAdd_

COPT = $(CC_LIB_FLAGS) $(COPTLESS)
C_OPT = $(CC_LIB_FLAGS) $(COPTLESS)
SPRNG_COMMON_DEPEND = primes_32.o primes_64.o \
			fwrap_mpi.o  cputime.o makeseed.o  store.o\
                    simple_mpi.o memory.o  communicate.o  checkid.o

OBJS = $(SPRNG_COMMON_DEPEND) sprng.o cmrg.o lcg.o lcg64.o lfg.o mlfg.o #pmlcg.o

all: $(TARGET)

libsprng.a: $(OBJS) 
	$(AR) $(ARFLAGS) libsprng.a $(OBJS) 
	$(RANLIB) libsprng.a

LIB_OBJS=$(OBJS:.o=_lo.o)

libsprng.$(LIB_SUFFIX): $(LIB_OBJS)
	$(DYNAMIC_LD) $(DYNAMIC_LD_FLAGS) libsprng.$(LIB_SUFFIX) $(LIB_FLAGS) $(LIB_OBJS) 

JUNK =  *~ core ii_files ti_files so_locations rii_files

clean:
	rm -fr $(LIB_OBJS) $(OBJS) $(TARGET) $(JUNK)

uninstall: 
	rm -fr $(LOCI_INSTALL_PATH)/lib/$(TARGET)

distclean:
	rm -fr $(OBJS) $(LIB_OBJS) $(TARGET) $(JUNK) #$(DEPEND_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif

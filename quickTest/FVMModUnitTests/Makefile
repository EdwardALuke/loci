# LOCI_BASE should be set before the Makefile
# Set TARGET to the name of your program
# Set FILES to list '.loci' files that will be compiled into your module, or

TARGET=FVMModUnitTest
FILES = $(AUTOMATIC_FILES)
INCLUDES = -I./include

###########################################################################
# No changes required below this line
###########################################################################

include $(LOCI_BASE)/Loci.conf

AUTOMATIC_FILES = $(call loci_compile_files,)
AUTOMATIC_OBJS = $(call loci_file2objs,$(FILES))

OBJS := $(AUTOMATIC_OBJS)


$(TARGET): $(OBJS)
	$(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

tests: $(TARGET)
	$(MAKE) -C tests TestResults LOCI_BASE=$(LOCI_BASE)

clean:
	rm -fr $(OBJS) $(TARGET)
	$(MAKE)  -C tests clean

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*)  core debug output $(wildcard *.o) 
# ".cc" files created from .loci files
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.$(LPP_I_SUFFIX))

distclean:
	rm -fr $(OBJS) $(TARGET) $(JUNK) $(LOCI_LPP_FILES) $(DEPEND_FILES)
	$(MAKE) -C tests distclean

DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
ifeq ($(filter $(MAKECMDGOALS),clean distclean ),)
-include $(DEPEND_FILES)
endif

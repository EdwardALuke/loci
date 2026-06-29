# Set up LOCI_BASE to point to the location where Loci is installed
LOCI_BASE ?= /simcenter/data4/loci/opt/Loci-Linux-i686-mpic++-rel-3-0-p8

# Include Loci.conf to get all of the rules for compiling Loci programs
include $(LOCI_BASE)/Loci.conf

TARGET = lfvm

# List your object files here
OBJS = main.o grid_reader.o setup_edges.o \
       gradient.o rhs.o integration.o interpolation.o stable.o \
       heat_bc.o edge_props.o cell_props.o physics.o area_reduce.o

# Junk files created by compilers that may need to be cleaned up
JUNK = *~ *.ii core ti_files 

default:
	  $(MAKE) $(TARGET)

# To compile a Loci program, use the .cc or .C prefix on your Loci source 
# modules and link using a command similar to below.
$(TARGET): $(OBJS) FRC
	  $(LD) -o $(TARGET) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC:

clean:
	  rm -fr $(OBJS) $(TARGET) $(JUNK) make.depend make.depend.bak


DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
include $(DEPEND_FILES)

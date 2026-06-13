# Set up LOCI_BASE to point to the location where Loci is installed
LOCI_BASE ?= /simcenter/data4/loci/opt/Loci-Linux-i686-mpic++-rel-3-0-p8

# Include Loci.conf to get all of the rules for compiling Loci programs
include $(LOCI_BASE)/Loci.conf

COPTLESS=-g
COPT=-g

TARGET = entities containers relations helpers

# List your object files here

EX1_OBJS = entities.o
EX2_OBJS = containers.o
EX3_OBJS = relations.o
EX4_OBJS = helpers.o

OBJS = $(EX1_OBJS) $(EX2_OBJS) $(EX3_OBJS) $(EX4_OBJS)

# Junk files created by compilers that may need to be cleaned up
JUNK = *~ *.ii core ti_files 

default:
	$(MAKE) $(TARGET)

# To compile a Loci program, use the .cc or .C prefix on your Loci source 
# modules and link using a command similar to below.
entities: $(EX1_OBJS)
	$(LD) -o entities $(EX1_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

containers: $(EX2_OBJS)
	$(LD) -o containers $(EX2_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

relations: $(EX3_OBJS)
	$(LD) -o relations $(EX3_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

helpers: $(EX4_OBJS)
	$(LD) -o helpers $(EX4_OBJS) $(LOCAL_LIBS) $(LIBS) $(LDFLAGS)

FRC:

clean:
	  rm -fr $(OBJS) $(TARGET) $(JUNK) make.depend make.depend.bak


DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
include $(DEPEND_FILES)

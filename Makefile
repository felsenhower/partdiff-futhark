# Futhark backend used to compile partdiff.fut: c, multicore, cuda, opencl, or hip.
FUTHARK_BACKEND ?= multicore

CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -O3 -g
LDFLAGS = $(CFLAGS)
LDLIBS  = -lm

ifeq ($(FUTHARK_BACKEND),multicore)
  LDLIBS += -lpthread
else ifeq ($(FUTHARK_BACKEND),opencl)
  LDLIBS += -lOpenCL
else ifeq ($(FUTHARK_BACKEND),cuda)
  LDLIBS += -lcuda -lcudart -lnvrtc
else ifeq ($(FUTHARK_BACKEND),hip)
  LDLIBS += -lhiprtc -lamdhip64
endif

all: partdiff

partdiff: partdiff.o partdiff_futhark.o

partdiff.o: partdiff.c partdiff_futhark.h

partdiff_futhark.c partdiff_futhark.h &: partdiff.fut
	futhark $(FUTHARK_BACKEND) --library -o partdiff_futhark partdiff.fut

clean:
	$(RM) partdiff.o partdiff_futhark.o partdiff_futhark.c partdiff_futhark.h partdiff_futhark.json partdiff

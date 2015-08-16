# Might fail to compile with -Wall
# CC = g++ -Wall

CC = g++

# CC = g++ -ansi -Wall -pedantic
# CFLAGS = -g -Wall -O3 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
# CFLAGS = -g -Wall

LDFLAGS = -lgsl -lm -lgslcblas

# these paths are not going to work for general audience
# although somehow the app compiles just fine with them
# still...
# GSL_INCLUDE = /home/chongw/include
# GSL_LIB = /home/chongw/lib

GSL_INCLUDE_MAC = /usr/local/include/
GSL_LIB_MAC = /usr/local/lib/

GSL_INCLUDE = GSL_INCLUDE_MAC
GSL_LIB = GSL_LIB_MAC

LSOURCE = main.cpp utils.cpp corpus.cpp ctr.cpp data.cpp opt.cpp
LHEADER = utils.h corpus.h ctr.h data.h opt.h

mac: $(LSOURCE) $(HEADER)
	  $(CC) -I$(GSL_INCLUDE_MAC) -L$(GSL_LIB_MAC) $(LSOURCE) -o ctr $(LDFLAGS)

mac-d: $(LSOURCE) $(HEADER)
	  $(CC) -g -I$(GSL_INCLUDE_MAC) -L$(GSL_LIB_MAC) $(LSOURCE) -o ctr $@ $(LDFLAGS)

linux: $(LSOURCE) $(HEADER)
	  $(CC) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(LSOURCE) -o ctr-condor $(LDFLAGS)

linux-d: $(LSOURCE) $(HEADER)
	  $(CC) -g -I$(GSL_INCLUDE) -L$(GSL_LIB_MAC) $(LSOURCE) -o ctr-condor $@ $(LDFLAGS)

clean:
	-rm -f ctr
clean-d:
	-rm -f ctr

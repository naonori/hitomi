
CC = g++
CFLAGS = -O3 -Wall
#CFLAGS += -D__MPI__

INCLUDE += -I/mwork0/sugiymnn/WORK/cosmo/fftw3/include/
INCLUDE += -I/mwork0/sugiymnn/WORK/cosmo/gsl/include/
INCLUDE += -I/mwork0/sugiymnn/WORK/cosmo/hitomi_measurement/

LIB += -L/mwork0/sugiymnn/WORK/cosmo/fftw3/lib/
LIB += -L/mwork0/sugiymnn/WORK/cosmo/gsl/lib/
LIB += -lm -lfftw3 -lgsl -lgslcblas
LIB += -std=c++11

PROGRAM = a.out

$(PROGRAM):main.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $< $(CLIBS) $(LIB) 

clean:	
	rm -f $(PROGRAM) core

distclean:	
	rm -f $(PROGRAM) *~


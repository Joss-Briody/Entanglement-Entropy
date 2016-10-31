# Define the source code modules for the library:

SOURCE = test1.cc

default: execs 

.PHONY: default execs output clean plots all

.SUFFIXES: .cc .hh .o .gp .cpp ,hpp

CC         = g++ -O2 -g -pedantic -Wall -I ./
INCLUDE    = -I$(CAV_INC) -I$(GSL_INC) -I$(FFTW_INC) 
LIBS       = $(CAV_LIB) $(GSL_LIB) $(FFTW_LIB) 
RM         = /bin/rm -f

%.o : %.cc
	$(CC) $(INCLUDE) -c -o $@ $<

%.o : %.cpp
	$(CC) $(INCLUDE) -c -o $@ $<

%.exe : %.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIBS)

%.exe : %.cpp
	$(CC) $(INCLUDE) -o $@ $< $(LIBS)

%: %.cc
	$(CC) $(INCLUDE) -o $@ $< $(LIBS)

%: %.cpp
	$(CC) $(INCLUDE) -o $@ $< $(LIBS)

%.out: %.exe
	./$<  $(DUMMYARGS) >  ./$@

execs: $(SOURCE:.cc=.exe) 

output: $(SOURCE:.cc=.out) 

clean:
	$(RM) $(SOURCE:.cc=.o) $(SOURCE:.cc=.exe) 

cleanout:
	$(RM) $(SOURCE:.cc=.out) 



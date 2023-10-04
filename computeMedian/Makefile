CPP = g++
CC = gcc
CFLAGS = -Wall  -O3  -fopenmp
LIBS = -lcfitsio -lm 
GLIBS = 
GLIBS += 
OBJECTS = computeMedian.o 
HEADERS = globalConstants.h

ALL : computeMedian.exe
	@echo "Listo!"

computeMedian.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o computeMedian.exe $(LIBS) $(GLIBS) $(CFLAGS)

computeMedian.o : computeMedian.cc $(HEADERS)
	$(CPP) -c computeMedian.cc -o computeMedian.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe

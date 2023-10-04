CPP = g++
CC = gcc
CFLAGS = $(shell root-config --cflags) -O3
LIBS = $(shell root-config --glibs) -lMinuit -lFoam
GLIBS =
GLIBS +=
OBJECTS = eventViewer.o listBox.o listBoxDict.o 
HEADERS = listBox.h 

ALL : eventViewer.exe
	echo "Listo!"

eventViewer.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o eventViewer.exe $(LIBS) $(GLIBS) $(CFLAGS)

eventViewer.o : eventViewer.cc $(HEADERS)
	$(CPP) -c eventViewer.cc -o eventViewer.o $(CFLAGS) $(LIBS)

listBox.o : listBox.cxx $(HEADERS)
	$(CPP) -c listBox.cxx -o listBox.o $(CFLAGS) $(LIBS)

listBoxDict.o : listBoxDict.C 
	$(CPP) -c listBoxDict.C -o listBoxDict.o $(CFLAGS) $(LIBS)

listBoxDict.C : listBox.h
	rootcint -f listBoxDict.C -c listBox.h

clean:
	rm -f *~ *.o *.exe listBoxDict.C listBoxDict.h 

CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -g -Wall -I$(CFITSIO) $(shell root-config --cflags) -O3
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --glibs) -lTreeViewer
GLIBS =
GLIBS +=
OBJECTS = main.o functions.o
HEADERS =
DEBUG_FLAG = -g

ALL : main.exe
	@echo File has been successfully compiled $(NEWLINE)

main.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) $(DEBUG_FLAG) -o main.exe $(LIBS) $(GLIBS) $(CFLAGS)

main.o : main.cpp $(HEADERS)
	$(CPP) -c main.cpp $(DEBUG_FLAG) -o main.o $(CFLAGS)

functions.o : functions.cpp $(HEADERS)
	$(CPP) -c functions.cpp $(DEBUG_FLAG) -o functions.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe

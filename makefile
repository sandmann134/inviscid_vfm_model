CC = g++
CFLAGS = -Wall -std=c++17 
VTK_LIBS = -lvtkCommonCore-9.2 -lvtkFiltersCore-9.2 -lvtkCommonDataModel-9.2 -lvtksys-9.2 -lvtkIOXML-9.2
INCLUDE_DIR = /usr/local/Cellar/vtk/9.2.6_5/include/vtk-9.2
LIB_DIR = /usr/local/Cellar/vtk/9.2.6_5/lib

all: test2

test2: test2.o filament.o
	$(CC) $(CFLAGS) -o test2 test2.o filament.o $(VTK_LIBS) -L $(LIB_DIR) -I $(INCLUDE_DIR)

test2.o: test2.cpp filament.h
	$(CC) $(CFLAGS) -c test2.cpp -I $(INCLUDE_DIR)

filament.o: filament.cpp filament.h
	$(CC) $(CFLAGS) -c filament.cpp -I $(INCLUDE_DIR)

clean:
	rm -f test2 test2.o filament.o

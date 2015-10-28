CC = gcc
CXX = g++
BOINC_DIR = ../../..
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_ZIP_DIR = $(BOINC_DIR)/zip
FREETYPE_DIR = /usr/include/freetype2
CFLAGS = -c -O3 $(ARCH) -Wall -Wextra -pedantic $(VARIANTFLAGS) -I$(BOINC_DIR) -I$(BOINC_LIB_DIR) -I$(BOINC_API_DIR) -I$(BOINC_ZIP_DIR) -I$(FREETYPE_DIR)
CXXFLAGS = $(CFLAGS) -std=c++11
LDFLAGS = $(ARCH) -L/usr/X11R6/lib -L.
CXXSOURCES = Graph.cpp boinc_functions.cpp utility.cpp pc.cpp
CSOURCES = erf.c
OBJECTS = $(CXXSOURCES:.cpp=.o) $(CSOURCES:.c=.o) BoincFile.o
EXECUTABLE = pc

all: $(CXXSOURCES) $(CSOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o ../bin/$@ -pthread libstdc++.a $(BOINC_API_DIR)/libboinc_api.a $(BOINC_LIB_DIR)/libboinc.a

BoincFile.o: BoincFile.cpp BoincFile.hpp
	$(CXX) $(CXXFLAGS) BoincFile.cpp -pthread -o BoincFile.o

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf ../bin/$(EXECUTABLE) *.o *~

.PHONY: all clean

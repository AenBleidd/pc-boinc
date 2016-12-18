CC = gcc
CXX = g++

BOINC_DIR = ../../..
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_ZIP_DIR = $(BOINC_DIR)/zip
BOINC_LIBS = $(BOINC_API_DIR)/libboinc_api.a $(BOINC_LIB_DIR)/libboinc.a

ifdef BOINC_STUB
BOINC_DIR = ../boinc_stub
BOINC_LIBS = 
endif

FREETYPE_DIR = /usr/include/freetype2
CFLAGS = -c -O3 $(ARCH) -Wall -Wextra -pedantic -Werror $(VARIANTFLAGS) -I$(BOINC_DIR) -I$(BOINC_LIB_DIR) -I$(BOINC_API_DIR) -I$(BOINC_ZIP_DIR) -I$(FREETYPE_DIR) -Isimd -MMD -MP
CXXFLAGS = $(CFLAGS) -std=c++11
LDFLAGS = $(ARCH) -L/usr/X11R6/lib -L.
CXXSOURCES = Graph.cpp boinc_functions.cpp utility.cpp pc.cpp
CSOURCES = erf.c
OBJECTS = $(CXXSOURCES:.cpp=.o) $(CSOURCES:.c=.o) BoincFile.o
EXECUTABLE = ../bin/pc

all: $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@ -static-libgcc -static-libstdc++ -static -pthread $(BOINC_LIBS)

BoincFile.o: BoincFile.cpp BoincFile.hpp
	$(CXX) $(CXXFLAGS) BoincFile.cpp -pthread -o BoincFile.o

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf ../bin/$(EXECUTABLE) *.o *~ *.d

.PHONY: all clean

-include $(CXXSOURCES:.cpp=.d) $(CSOURCES:.c=.d) BoincFile.d

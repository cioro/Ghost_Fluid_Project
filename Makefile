CC=g++
CFLAGS=-c -Wall -Wextra -std=c++11 -g -Wshadow
LDFLAGS= 
LIBS=  -lm

SOURCES:=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE = output

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)
.cpp.o:
	$(CC) $(CFLAGS) $(LIBS) $< -o $@
clean:
		rm -rf *.o

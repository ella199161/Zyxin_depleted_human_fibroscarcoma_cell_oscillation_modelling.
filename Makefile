#.SUFFIXES:.cpp .o
CC=g++
CFLAGS= -Wall
CFLAGS+= -O3 
CFLAGS+= -std=c++11

#CFLAGS+= -gstabs+
#CFLAGS+= -pg


LINKER= $(CC)
LFLAGS= $(CFLAGS)

SRC = set1.cpp Source.cpp

OBJ = $(SRC:.cpp=.o)

HEADERS= parameter.h Source.h

EXEC = set1

%.o:	%.cpp $(HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

all:	$(HEADERS) $(EXEC)

$(EXEC): 	$(OBJ)
		$(LINKER) $(LFLAGS) $(OBJ) -o $@

clean:
	rm -f *.o; rm -f $(EXEC); rm -f core; rm -f *.out; rm -f *~;




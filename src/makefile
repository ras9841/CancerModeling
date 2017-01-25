OBJS = main.o cells.o simulation.o
CC = g++
CDEBUG = -g
CFLAGS = -Wall -std=c++11 -c $(CDEBUG)
LFLAGS = -Wall -lconfig++
DEBUG = off

run : $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o run

main.o : main.cpp params.hpp 
ifeq ($(DEBUG), off)
	$(CC) $(CFLAGS) main.cpp params.hpp
else
	$(CC) $(CFLAGS) -DDEBUG main.cpp params.hpp
endif

cells.o : cells.cpp cells.hpp 
ifeq ($(DEBUG), off)
	$(CC) $(CFLAGS) cells.cpp cells.hpp
else
	$(CC) $(CFLAGS) -DDEBUG cells.cpp cells.hpp
endif

simulation.o : simulation.cpp simulation.hpp 
ifeq ($(DEBUG), off)
	$(CC) $(CFLAGS) simulation.cpp simulation.hpp
else
	$(CC) $(CFLAGS) -DDEBUG simulation.cpp simulation.hpp
endif

clean:
	rm *.o *.gch run 

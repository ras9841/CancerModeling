OBJS = main.o cells.o simulation.o table.o
CC = g++
CDEBUG = -g
CFLAGS = -Wall -std=c++11 -c $(CDEBUG)
LFLAGS = -Wall -lconfig++
DEBUG = off
CELL = cells.cpp cells.hpp
TABLE = table.cpp table.hpp
SIM = simulation.cpp simulation.hpp

run : $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o run

opt : $(OBJS) 
	$(CC) -O3 $(LFLAGS) $(OBJS) -o run

main.o : main.cpp params.hpp $(SIM)
ifeq ($(DEBUG), off)
	$(CC) -O3 $(CFLAGS) main.cpp params.hpp $(SIM)
else
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" main.cpp params.hpp $(SIM)
endif

cells.o : $(CELL)
ifeq ($(DEBUG), off)
	$(CC) -O3 $(CFLAGS) $(CELL)
else
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" $(CELL)
endif

simulation.o : $(CELL) $(SIM) $(TABLE)
ifeq ($(DEBUG), off)
	$(CC) -O3 $(CFLAGS) $(CELL) $(SIM) $(TABLE)
else
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" $(CELL) $(SIM) $(TABLE)
endif

table.o : $(CELL) $(SIM) $(TABLE)
ifeq ($(DEBUG), off)
	$(CC) -O3 $(CFLAGS) $(CELL) $(SIM) $(TABLE)
else
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" $(CELL) $(SIM) $(TABLE)
endif

clean:
	rm *.o *.gch run 

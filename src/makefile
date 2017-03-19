OBJS = main.o cells.o simulation.o table.o
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
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" main.cpp params.hpp
endif

cells.o : cells.cpp cells.hpp 
ifeq ($(DEBUG), off)
	$(CC) $(CFLAGS) cells.cpp cells.hpp
else
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" cells.cpp cells.hpp
endif

simulation.o : simulation.cpp simulation.hpp table.hpp
ifeq ($(DEBUG), off)
	$(CC) $(CFLAGS) simulation.cpp simulation.hpp
else
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" simulation.cpp simulation.hpp
endif

table.o : table.cpp table.hpp linked_list.hpp simulation.hpp
ifeq ($(DEBUG), off)
	$(CC) $(CFLAGS) table.cpp table.hpp linked_list.hpp simulation.hpp
else
	$(CC) $(CFLAGS) -DDEBUG=\"$(DEBUG)\" table.cpp table.hpp linked_list.hpp simulation.hpp

endif

clean:
	rm *.o *.gch run 

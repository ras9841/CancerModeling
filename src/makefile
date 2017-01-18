OBJS = main.o
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

clean:
	rm *.o *.gch run 

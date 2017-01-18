OBJS = main.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -std=c++11 -c $(DEBUG)
LFLAGS = -Wall -lconfig++

run : $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o run

main.o : main.cpp params.hpp 
	$(CC) $(CFLAGS) main.cpp params.hpp

clean:
	rm *.o *.gch run

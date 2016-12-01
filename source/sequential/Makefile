CC=g++
FLAGS=-lgmpxx -lgmp
EXEC=main

ALL: factors.o batchGCD.o main.o
	$(CC) -o $(EXEC) factors.o batchGCD.o main.o $(FLAGS)

main.o: main.cpp factors.cpp factors.hpp batchGCD.cpp batchGCD.hpp
	$(CC) -c main.cpp factors.cpp batchGCD.cpp

factors.o: factors.cpp factors.hpp
	$(CC) -c factors.cpp

batchGCD.o: batchGCD.cpp batchGCD.hpp
	$(CC) -c batchGCD.cpp

clean:
	rm -rf test
	rm -f *~
	rm -rf *.o

test_tiny:
	./$(EXEC) 100 ../../data/tiny.dat

test_big:
	./$(EXEC) 2000 ../../data/gen512-3233331101-0.dat

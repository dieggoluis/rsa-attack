CC=mpic++
FLAGS=-lgmpxx -lgmp
EXEC=main
NPROC_TINY=20
NPROC_BIG=50

ALL: factors.o batchGCD.o main.o
	$(CC) -o $(EXEC) factors.o batchGCD.o main.o $(FLAGS)

main.o: main.cpp factors.cpp factors.hpp batchGCD.cpp batchGCD.hpp
	$(CC) -c main.cpp factors.cpp batchGCD.cpp

factors.o: factors.cpp factors.hpp
	$(CC) -c factors.cpp

batchGCD.o: batchGCD.cpp batchGCD.hpp
	$(CC) -c batchGCD.cpp

clean:
	rm -rf main
	rm -rf *.o
	rm -f *~

test_tiny:
	mpirun -np $(NPROC_TINY) $(EXEC) ../../data/tiny.dat

test_big:
	mpirun -np $(NPROC_BIG) $(EXEC) ../../data/gen512-3233331101-0.dat

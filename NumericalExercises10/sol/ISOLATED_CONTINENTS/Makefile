CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

tsp.exe : tsp.o specimen.o population.o vertices.o random.o
	$(CC) specimen.o population.o vertices.o tsp.o random.o -o tsp.exe
tsp.o : tsp.cpp tsp.h
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
vertices.o : vertices.cpp 
	$(CC) -c vertices.cpp -o vertices.o $(CFLAGS)
population.o : population.cpp specimen.cpp 
	$(CC) -c population.cpp -o population.o
specimen.o: specimen.cpp 
	$(CC) -c specimen.cpp -o specimen.o
random.o : random.cpp
	$(CC) -c random.cpp -o random.o
clean :
	rm *.o tsp.exe seed.out
run :
	mpiexec -np 8 ./tsp.exe

CC = g++
CFLAGS = -Wall -O3 --std=c++11

SAN.exe : SAN.o random.o
	$(CC) random.o SAN.o -o SAN.exe
SAN.o : SAN.cpp SAN.h
	$(CC) -c SAN.cpp -o SAN.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SAN.exe seed.out

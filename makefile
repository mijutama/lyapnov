
all: run

clean:; rm run main.o

run: main.o
	g++ main.o -o run -lm -O2

main.o: main.cpp
	g++ -c main.cpp -o main.o -std=c++11


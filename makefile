
all: run fft reat branch

clean:; rm run main.o

run: main.o
	g++ main.o -o run -lm -O2

main.o: main.cpp
	g++ -c main.cpp -o main.o -std=c++11

fft: fft.cpp
	g++ fft.cpp -o fft -O2 -std=c++11 -lfftw3

reat: reat.cpp
	g++ reat.cpp -o reat -O2 -std=c++11

branch: branch.cpp
	g++ branch.cpp -o branch -O2 -std=c++11


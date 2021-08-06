all: filterChimericBam

filterChimericBam.o: filterChimericBam.cpp
	g++ -Wall -std=c++11 -c filterChimericBam.cpp -o filterChimericBam.o

filterChimericBam: filterChimericBam.o
	g++ filterChimericBam.o -o filterChimericBam

clean:
	rm *.o filterChimericBam
all: filterChimericBam getCommonReads sortLikeSamtools

filterChimericBam.o: filterChimericBam.cpp
	g++ -Wall -std=c++11 -c filterChimericBam.cpp -o filterChimericBam.o

filterChimericBam: filterChimericBam.o
	g++ filterChimericBam.o -o filterChimericBam

getCommonReads.o: getCommonReads.cpp
	g++ -Wall -std=c++11 -c getCommonReads.cpp -o getCommonReads.o

getCommonReads: getCommonReads.o
	g++ getCommonReads.o -o getCommonReads

sortLikeSamtools.o: sortLikeSamtools.cpp
	g++ -Wall -std=c++11 -c sortLikeSamtools.cpp -o sortLikeSamtools.o

sortLikeSamtools: sortLikeSamtools.o
	g++ sortLikeSamtools.o -o sortLikeSamtools

clean:
	rm *.o filterChimericBam getCommonReads sortLikeSamtools
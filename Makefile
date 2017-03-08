all: code-mine.cpp
	g++ -std=c++11 -Ofast -g -o mypart code-mine.cpp

clean:
	rm -f mypart*.o

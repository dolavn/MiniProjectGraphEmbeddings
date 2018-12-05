all: graphEmbed

graphEmbed: bin/Main.o bin/Graph.o bin/GraphGenerator.o
	g++ -o bin/graphEmbed bin/Main.o bin/Graph.o bin/GraphGenerator.o
	
bin/Main.o: src/Main.cpp
	g++ -g -Wall -Weffc++ -std=c++11 -c -Iinclude -o bin/Main.o src/Main.cpp

bin/Graph.o: src/Graph.cpp
	g++ -g -Wall -Weffc++ -std=c++11 -c -Iinclude -o bin/Graph.o src/Graph.cpp

bin/GraphGenerator.o: src/GraphGenerator.cpp
	g++ -g -Wall -Weffc++ -std=c++11 -c -Iinclude -o bin/GraphGenerator.o src/GraphGenerator.cpp


clean:
	rm -f bin/*
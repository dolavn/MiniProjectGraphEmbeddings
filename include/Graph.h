#include <vector>
#include <string>

#ifndef GRAPH_H_
#define GRAPH_H_

class Node;
class Edge;

typedef struct indDistPair{
	int nodeInd;
	double dist;
	
	indDistPair(int nodeInd,double dist):nodeInd(nodeInd),dist(dist){}
	indDistPair():nodeInd(0),dist(0){}
} indDistPair;

typedef struct preEdge{
    int ind1;
    int ind2;
    double weight;
    
    preEdge(int ind1,int ind2,double weight):ind1(ind1),ind2(ind2),weight(weight){}
} preEdge;

typedef struct triangle{
    int ind1;
    int ind2;
    int ind3;
    triangle(int ind1,int ind2,int ind3):ind1(ind1),ind2(ind2),ind3(ind3){}
} triangle;

bool cmp(Edge* edge1, Edge* edge2);
bool cmpHeap(indDistPair a, indDistPair b);

bool freePairingExists(std::vector<int>* pairing,unsigned int deg,unsigned int ind);

bool goodPairing(std::vector<int>* pairing,int i,int j,int deg);

std::vector<int>* createPairing(unsigned int size,unsigned int deg);

std::string printPair(indDistPair a);

preEdge translateLine(std::string tString);

class Graph {
public:
	static Graph* getRandomGraph(int size, double p, double minWeight,double maxWeight);
	static Graph* getBiPartiteGraph(int size1,int size2, double p,double minWeight,double maxWeight);
	static Graph* getRegularGraph(int size, int deg, double minWeight,double maxWeight);
	static Graph* getPlanarGraph(int size, double p,double minWeight,double maxWeight);
        static std::vector<Graph*> readFromFile(std::string path);
	static bool checkIfSpanner(Graph& g, Graph& spanner,double t);
	Graph();
	Graph(int n);
	~Graph();
	Graph* createSpanner(double r);
        Graph* createInverse(double minWeight,double maxWeight);
        Graph* createMST();
	std::vector<double> shortestPathVec(int i, int ind2);
	std::vector<double> shortestPath(int i);
        std::vector<Edge*> getEdges(){return edges;}
	double shortestPath(int ind1, int ind2);
        bool onSameComponent(int ind1,int ind2);
	double getGraphWeight();
	void addNode();
        int numOfNodes(){return nodes.size();}
	float getAverageDegree(){return ((float)edges.size()*2)/(float)nodes.size();}
	void printMatrix();
	void printEdges();
	void addEdge(int ind1, int ind2, double weight);
        void removeEdge(int ind1,int ind2);
	void sortEdges();
	void printDistances();
	Node& getNode(int ind);
	int numOfEdges();
        void saveToFile(std::ofstream& fileStream);
private:
        static Graph* getBasePlanarGraph(int size,double minWeight,double maxWeight);
        void dfsVisit(int ind,std::vector<bool>* visited);
	std::vector<Node*> nodes;
	std::vector<Edge*> edges;
	std::vector<std::vector<Edge*>> matrix;
};

class Node {
public:
	Node(Graph& graph, int ind);
	void addNeighbour(int ind, double weight);
	void addEdge(Edge* edge);
        void removeEdge(int otherInd);
	std::vector<Edge*> getNeighbours();
	double shortestPathTo(int other);
        bool onSameComponent(int other);
        unsigned int getDegree();
private:
	Graph& graph;
	int ind;
	std::vector<Edge*> neighbours;
};

class Edge {
public:
	Edge(Graph& graph, int ind1, int ind2, double weight);
	std::string toString();
	double getWeight();
	int getOther(int i);
	int getInd1();
	int getInd2();
private:
	Graph& graph;
	int ind1;
	int ind2;
	double weight;
};

#endif
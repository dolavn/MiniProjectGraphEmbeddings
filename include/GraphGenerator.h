#ifndef GRAPH_GENERATOR_H_
#define GRAPH_GENERATOR_H_
#include <string>
class Graph;

class GraphGenerator{
public:
    virtual Graph* generateGraph()=0;
    virtual std::string toString()=0;
    virtual ~GraphGenerator();
};

class RandomGraphGenerator:public GraphGenerator{
public:
    RandomGraphGenerator(int size,double p,double minWeight,double maxWeight);
    Graph* generateGraph();
    std::string toString();
    virtual ~RandomGraphGenerator();
private:
    int size;
    double p;
    double minWeight;
    double maxWeight;
};

class RegularGraphGenerator:public GraphGenerator{
public:
    RegularGraphGenerator(int size,int deg,double minWeight,double maxWeight);
    Graph* generateGraph();
    std::string toString();
    virtual ~RegularGraphGenerator();
private:
    int size;
    int deg;
    double minWeight;
    double maxWeight;
};

class BipartiteGraphGenerator:public GraphGenerator{
public:
    BipartiteGraphGenerator(int size1,int size2,double p,double minWeight,double maxWeight);
    Graph* generateGraph();
    std::string toString();
    virtual ~BipartiteGraphGenerator();
private:
    int size1;
    int size2;
    double p;
    double minWeight;
    double maxWeight;
};

class PlanarGraphGenerator:public GraphGenerator{
public:
    PlanarGraphGenerator(int size,double p,double minWeight,double maxWeight);
    Graph* generateGraph();
    std::string toString();
    virtual ~PlanarGraphGenerator(); 
private:
    int size;
    double p;
    double minWeight;
    double maxWeight;
};

class FileGraphGenerator:public GraphGenerator{
public:
    FileGraphGenerator(std::string path);
    Graph* generateGraph();
    std::string toString();
    virtual ~FileGraphGenerator();
    int size(){return graphs.size();}
private:
    std::string path;
    std::vector<Graph*> graphs;
    unsigned int curr;
};


#endif
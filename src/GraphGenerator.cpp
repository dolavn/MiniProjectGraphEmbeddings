#include "Graph.h"
#include "GraphGenerator.h"
#include <stdexcept>
#include <sstream>

GraphGenerator::~GraphGenerator(){
}

RandomGraphGenerator::RandomGraphGenerator(int size,double p,double minWeight,double maxWeight):size(size),p(p),minWeight(minWeight),maxWeight(maxWeight){
    if(p<0 || p>1){
        throw std::runtime_error("Illegel P value");
    }
}

Graph* RandomGraphGenerator::generateGraph(){
    return Graph::getRandomGraph(size,p,minWeight,maxWeight);
}

std::string RandomGraphGenerator::toString(){
    std::stringstream ss;
    ss << "random " << size << " nodes graphs" << " p=" << p << " weight function between " << minWeight <<" and " << maxWeight;
    return ss.str();
}

RandomGraphGenerator::~RandomGraphGenerator(){
}

RegularGraphGenerator::RegularGraphGenerator(int size,int deg,double minWeight,double maxWeight):size(size),deg(deg),minWeight(minWeight),maxWeight(maxWeight){
    if(size*deg%2!=0){
        throw std::runtime_error("Sum of degrees is uneven. Unable to construct");
    }
    if(deg>=size){
        throw std::runtime_error("Degree of each node is larger than number of nodes. Unable to construct");
    }
}

Graph* RegularGraphGenerator::generateGraph(){
    return Graph::getRegularGraph(size,deg,minWeight,maxWeight);
}

std::string RegularGraphGenerator::toString(){
    std::stringstream ss;
    ss << deg << "-regular " << size << " nodes graphs weight function between " << minWeight << " and " << maxWeight;
    return ss.str();
}

RegularGraphGenerator::~RegularGraphGenerator(){
}

BipartiteGraphGenerator::BipartiteGraphGenerator(int size1,int size2,double p,double minWeight,double maxWeight):size1(size1),size2(size2),p(p),minWeight(minWeight),maxWeight(maxWeight){
    if(p<0 || p>1){
        throw std::runtime_error("Illegel P value");
    }
}

Graph* BipartiteGraphGenerator::generateGraph(){
    return Graph::getBiPartiteGraph(size1,size2,p,minWeight,maxWeight);
}

std::string BipartiteGraphGenerator::toString(){
    std::stringstream ss;
    ss << "random bipartite graphs |V1|=" << size1 << " |V2|=" << size2 << " p=" << p << " weight function between " << minWeight << " and " << maxWeight;
    return ss.str();
}

BipartiteGraphGenerator::~BipartiteGraphGenerator(){
}

PlanarGraphGenerator::PlanarGraphGenerator(int size,double p,double minWeight,double maxWeight):size(size),p(p),minWeight(minWeight),maxWeight(maxWeight){
    if(p<0 || p>1){
        throw std::runtime_error("Illegel P value");
    }
}

Graph* PlanarGraphGenerator::generateGraph(){
    return Graph::getPlanarGraph(size,p,minWeight,maxWeight);
}

std::string PlanarGraphGenerator::toString(){
    std::stringstream ss;
    ss << "random " << size << " nodes planar graph, P=" << p << " weight function between " << minWeight << " and " << maxWeight;
    return ss.str();
}

PlanarGraphGenerator::~PlanarGraphGenerator(){
}

FileGraphGenerator::FileGraphGenerator(std::string path):path(path),graphs(),curr(0){
    graphs = Graph::readFromFile(path);
}

Graph* FileGraphGenerator::generateGraph(){
    return graphs[curr++];
}

std::string FileGraphGenerator::toString(){
    std::stringstream ss;
    ss << "graphs from file " << path;
    return ss.str();
}

FileGraphGenerator::~FileGraphGenerator(){
    for(unsigned int i=0;i<graphs.size();i++){
        delete(graphs[i]);
    }
}
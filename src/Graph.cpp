#include "Graph.h"
#include "Heap.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <time.h>
using namespace std;

Graph::Graph() :nodes(), edges(), matrix() {
}

Graph* Graph::getRandomGraph(int size, double p, double minWeight,double maxWeight) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
	srand((time_t)ts.tv_nsec);
	Graph* ans = new Graph(size);
	for (int i = 0; i<size; i++) {
		for (int j = i + 1; j<size; j++) {
			double r = ((double)rand() / (RAND_MAX));
			if (p>r) {
				double weight = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
				ans->addEdge(i, j, weight);
			}
		}
	}
	return ans;
}

Graph* Graph::getBiPartiteGraph(int size1,int size2,double p,double minWeight,double maxWeight){
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
	srand((time_t)ts.tv_nsec);
	Graph* ans = new Graph(size1+size2);
	for(int i=0;i<size1;i++){
		for(int j=size1;j<size1+size2;j++){
				double r = ((double)rand()/(RAND_MAX));
				if(p>r){
					double weight = (((int)rand())%((int)(maxWeight-minWeight+1)))+minWeight;
					ans->addEdge(i,j,weight);
				}
		}
	}
	return ans;
}

Graph* Graph::getRegularGraph(int size, int deg,double minWeight,double maxWeight) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
	srand((time_t)ts.tv_nsec);
        if(size*deg%2!=0 || deg>=size){
            throw 2;
        }
	Graph* ans = new Graph(size);
        bool inverse=false;
        if(deg>size/2){
            inverse=true;
            deg = size-1-deg;
        }
        vector<int>* points = createPairing(size,deg);
	for (int i = 0; i < size; i++) {
		for (int j = i*deg; j < (i+1)*deg; j++) {
			int node = (*points)[j]/deg; 
                        if(node>i){
                            double weight = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
                            ans->addEdge(i,node,weight);
                        }
		}
	}
	delete(points);
        if(inverse){
            Graph* inverse = ans->createInverse(minWeight,maxWeight);
            delete(ans);
            ans = inverse;
        }
	return ans;
}

Graph* Graph::createInverse(double minWeight,double maxWeight){
    Graph* ans = getRandomGraph(nodes.size(),1,minWeight,maxWeight); /*Creates full graph*/
    for(unsigned int i=0;i<nodes.size();i++){
        Node& currNode = getNode(i);
        for(unsigned int j=0;j<currNode.getNeighbours().size();j++){
            int otherInd = currNode.getNeighbours()[j]->getOther(i);
            ans->removeEdge(i,otherInd);
        }
    }
    return ans;
}

Graph* Graph::getPlanarGraph(int size,double p,double minWeight,double maxWeight) {
    Graph* ans = getBasePlanarGraph(size,minWeight,maxWeight);
    std::vector<Edge*> edges = ans->getEdges();
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    srand((time_t)ts.tv_nsec);
    for(unsigned int i=0;i<edges.size();i++){
        double r = ((double)rand() / (RAND_MAX));
        if(p<r){
            int ind1 = edges[i]->getInd1();
            int ind2 = edges[i]->getInd2();
            double weight = edges[i]->getWeight();
            ans->removeEdge(ind1,ind2);
            if(!ans->onSameComponent(ind1,ind2)){
                ans->addEdge(ind1,ind2,weight);
            }
        }
    }
    return ans;
}

Graph* Graph::getBasePlanarGraph(int size,double minWeight,double maxWeight){
    if(size<=3){
        return getRandomGraph(size,1,minWeight,maxWeight);
    }
    Graph* ans = new Graph(size);
    std::vector<triangle> faces = std::vector<triangle>();
    double weight1 = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
    double weight2 = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
    double weight3 = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
    ans->addEdge(0,1,weight1);
    ans->addEdge(1,2,weight2);
    ans->addEdge(0,2,weight3);
    triangle face1(0,1,2);
    faces.push_back(face1);
    for(int i=3;i<size;i++){
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        srand((time_t)ts.tv_nsec);
        int faceInd = ((int)(rand()))%faces.size();
        triangle currFace = faces[faceInd];
        weight1 = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
        weight2 = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
        weight3 = (((int)rand()) % ((int)(maxWeight-minWeight+1))) + minWeight;
        ans->addEdge(i,currFace.ind1,weight1);
        ans->addEdge(i,currFace.ind2,weight2);
        ans->addEdge(i,currFace.ind3,weight3);
        triangle newFace1(i,currFace.ind1,currFace.ind2);
        triangle newFace2(i,currFace.ind2,currFace.ind3);
        triangle newFace3(i,currFace.ind3,currFace.ind1);
        faces.push_back(newFace1);faces.push_back(newFace2);faces.push_back(newFace3);
        faces[faceInd] = faces[faces.size()-1];
        faces.pop_back();
    }
    return ans;
}

bool Graph::checkIfSpanner(Graph& g, Graph& spanner,double t) {
	if (g.nodes.size() != spanner.nodes.size()) { return false; }
	for (unsigned int i = 0; i < g.nodes.size(); i++) {
            for (unsigned int j = i + 1; j < g.nodes.size(); j++) {
                double distOrig = g.nodes[i]->shortestPathTo(j);
                double distNew = spanner.nodes[i]->shortestPathTo(j);
                if (distOrig*t < distNew) { 
                        return false;
                }
            }
	}
	return true;
}

Graph::~Graph() {
	for (unsigned int i = 0; i<nodes.size(); i++) {
		delete(nodes[i]);
	}
	for (unsigned int i = 0; i<edges.size(); i++) {
		delete(edges[i]);
	}
}

Graph::Graph(int n) :nodes(), edges(), matrix() {
	for (int i = 0; i<n; i++) {
		Node* newNode = new Node(*this, i);
		nodes.push_back(newNode);
		vector<Edge*> newCol = vector<Edge*>(n);
		matrix.push_back(newCol);
	}
}

void Graph::addNode() {
	Node* newNode = new Node(*this, nodes.size());
	nodes.push_back(newNode);
	for (unsigned int i = 0; i<matrix.size(); i = i + 1) {
		matrix[i].push_back(nullptr);
	}
	vector<Edge*> newCol = vector<Edge*>(nodes.size());
	matrix.push_back(newCol);
}

Graph* Graph::createSpanner(double r) {
	Graph* ans = new Graph(nodes.size()); //G'
	sortEdges();
	for (unsigned int i = 0; i<edges.size(); i++) { //for every edge [u,v] in G   
		Edge& currEdge = *(edges[i]);
		int ind1 = currEdge.getInd1();
		int ind2 = currEdge.getInd2();
		double currDist = ans->shortestPath(ind1, ind2);
		double weight = currEdge.getWeight();
		if (weight*r<currDist) {
                    ans->addEdge(ind1, ind2, weight);
		}
	}
	return ans;
}

Graph* Graph::createMST(){
    Graph* ans = new Graph(nodes.size());
    sortEdges();
    for(unsigned int i=0;i<edges.size();i++){
        Edge& currEdge = *(edges[i]);
        int ind1 = currEdge.getInd1();
        int ind2 = currEdge.getInd2();
        double weight = currEdge.getWeight();
        if(!(ans->onSameComponent(ind1,ind2))){
            ans->addEdge(ind1,ind2,weight);
        }
    }
    return ans;
}

void Graph::saveToFile(std::ofstream& fileStream){
    fileStream << "G\n"; //signifies new graph
    fileStream << nodes.size() << "\n";
    for(unsigned int i=0;i<edges.size();i++){
        Edge& currEdge = *edges[i];
        int ind1 = currEdge.getInd1();
        int ind2 = currEdge.getInd2();
        double weight = currEdge.getWeight();
        fileStream << ind1 << "," << ind2 << "," << weight << "\n";
    }
    fileStream << "E\n"; //signifes end of graph
}

std::vector<Graph*> Graph::readFromFile(std::string filePath){
    std::vector<Graph*> ans = std::vector<Graph*>();
    std::ifstream stream;
    std::string line;
    stream.open(filePath);
    if(stream.is_open()){
        while(getline(stream,line)){
            if(line=="G"){
                getline(stream,line);
                int size = std::stoi(line);
                Graph* curr = new Graph(size);
                while(getline(stream,line) && line!="E"){
                    preEdge currEdge = translateLine(line);
                    curr->addEdge(currEdge.ind1,currEdge.ind2,currEdge.weight);
                }
                ans.push_back(curr);
            }
        }
    }else{
        throw std::runtime_error("Couldn't open file!");
    }
    stream.close();
    return ans;
}

int getPInd(int pInd, int ind) {
	if (pInd == ind) { return 0; }
	if (pInd < ind) { return pInd + 1; }
	else{ return pInd; }
}

vector<double> Graph::shortestPathVec(int ind, int target) {
	vector<double> ans = vector<double>(nodes.size());\
	for (unsigned int i = 0; i<ans.size(); i++) {
		ans[i] = numeric_limits<double>::infinity(); //ans[i]=inf
	}
	ans[ind] = 0; //ans[source]=0
	vector<indDistPair> heapTemplate;
	heapTemplate.push_back(indDistPair(ind, 0));
	for (unsigned int i = 0; i<ans.size(); i++) {
		if (i != (unsigned int)ind) {
			indDistPair curr(i, ans[i]);
			heapTemplate.push_back(curr);
		}
	}
	Heap<indDistPair> heap = Heap<indDistPair>(heapTemplate, cmpHeap, printPair,true);
	vector<int*> pointers = heap.getPointerVec();
	int count = 0;
	while (heap.getSize()>0) { //dijkstra main loop
		count++;
		indDistPair removed = heap.popHead();
		int indMin = removed.nodeInd;
		if (indMin == target) { break; }
		int pInd = getPInd(indMin, ind);
		pointers[pInd] = nullptr;  //memory management
		Node& node = getNode(indMin);
		vector<Edge*> neighbours = node.getNeighbours();
		for (unsigned int i = 0; i<neighbours.size(); i++) {
			Edge& currEdge = *(neighbours[i]);
			int otherInd = currEdge.getOther(indMin); //The neighbour
			double newDist = ans[indMin] + currEdge.getWeight();
			if (newDist<ans[otherInd]) {
				ans[otherInd] = newDist;
				int otherHeapInd = *pointers[getPInd(otherInd,ind)];
				indDistPair newDistPair(otherInd,newDist); 
				heap.setValue(otherHeapInd, newDistPair);
			}
		}
	}
	return ans;
}

vector<double> Graph::shortestPath(int ind) {
	return shortestPathVec(ind, -1);
}

double Graph::shortestPath(int ind1, int ind2) {
	vector<double> distances = shortestPathVec(ind1, ind2);
	return distances[ind2];
}

bool Graph::onSameComponent(int ind1,int ind2){
    std::vector<bool> visited = std::vector<bool>(nodes.size());
    dfsVisit(ind1,&visited);
    return visited[ind2];
}


void Graph::dfsVisit(int ind1,std::vector<bool>* visited){
    Node& curr = getNode(ind1);
    (*visited)[ind1]=true;
    std::vector<Edge*> neighbours = curr.getNeighbours();
    for(unsigned int i=0;i<neighbours.size();i++){
        Edge& currEdge = *(neighbours[i]);
        int otherInd = currEdge.getOther(ind1);
        if(!(*visited)[otherInd]){
            dfsVisit(otherInd,visited);
        }
    }
}

double Graph::getGraphWeight() {
	double ans = 0;
	for (unsigned int i = 0; i<edges.size(); i++) {
		ans = ans + edges[i]->getWeight();
	}
	return ans;
}

void Graph::printMatrix() {
	for (unsigned int i = 0; i < matrix.size(); i++) {
		vector<Edge*>& col = matrix[i];
		for (unsigned j = 0; j < col.size(); j++) {
			string cell;
			if (nullptr == col[j]) {
				cell = "|  0   |";
			}
			else {
				cell = "|" + col[j]->toString() + "|";
			}
			cout << cell + " ";
		}
		cout << endl;
	}
}

void Graph::printEdges() {
	for (unsigned int i = 0; i<edges.size(); i++) {
		cout << edges[i]->toString() << "\n";
	}
	cout << endl;
}

void Graph::printDistances() {
	for (unsigned int i = 0; i < nodes.size(); i++) {
		for (unsigned int j = i + 1; j < nodes.size(); j++) {
			double dist = nodes[i]->shortestPathTo(j);
			std::cout << i << " to " << j << " " << dist << std::endl;
		}
	}
}

void Graph::sortEdges() {
    sort(edges.begin(), edges.end(), cmp);
}

void Graph::addEdge(int ind1, int ind2, double weight) {
	Edge* newEdge = new Edge(*this, ind1, ind2, weight);
	edges.push_back(newEdge);
	matrix[ind1][ind2] = newEdge;
	matrix[ind2][ind1] = newEdge;
	nodes[ind1]->addEdge(newEdge); nodes[ind2]->addEdge(newEdge);
}

void Graph::removeEdge(int ind1,int ind2){
    Edge* edge = matrix[ind1][ind2];
    if(edge==nullptr){
        return;
    }
    int delInd=-1;
    for(unsigned int i=0;i<edges.size() && delInd==-1;i++){
        if(edges[i]==edge){
            delInd=i;
        }
    }
    edges[delInd] = edges[edges.size()-1];
    edges.pop_back();
    Node& node1 = getNode(ind1);
    Node& node2 = getNode(ind2);
    matrix[ind1][ind2] = nullptr;
    matrix[ind2][ind1] = nullptr;
    node1.removeEdge(ind2);
    node2.removeEdge(ind1);
    delete(edge);
}

Node& Graph::getNode(int ind) {
	return *(nodes[ind]);
}

int Graph::numOfEdges() { return edges.size(); }

Node::Node(Graph& graph, int ind) :graph(graph), ind(ind), neighbours() {
}

void Node::addNeighbour(int ind, double weight) {
	graph.addEdge(this->ind, ind, weight);
}

void Node::removeEdge(int otherInd){
    for(unsigned int i=0;i<neighbours.size();i++){
        Edge* curr = neighbours[i];
        if(curr->getOther(ind)==otherInd){
            neighbours[i] = neighbours[neighbours.size()-1];
            neighbours[neighbours.size()-1] = curr;
            neighbours.pop_back();
            return;
        }
    }
}

void Node::addEdge(Edge* edge) {
	neighbours.push_back(edge);
}

vector<Edge*> Node::getNeighbours() {
	return neighbours;
}

bool Node::onSameComponent(int other){
    return graph.onSameComponent(ind,other);
}

double Node::shortestPathTo(int other) {
	return graph.shortestPath(ind, other);
}

unsigned int Node::getDegree(){
    return neighbours.size();
}

Edge::Edge(Graph& graph, int ind1, int ind2, double weight) :graph(graph), ind1(ind1), ind2(ind2), weight(weight) {
}

int Edge::getOther(int i) {
	if (i == ind1) {
		return ind2;
	}
	if (i == ind2) {
		return ind1;
	}
	return -1;
}

string Edge::toString() {
	stringstream ss;
	ss << "(" << ind1 << "," << ind2 << ")" << weight;
	return ss.str();
}

double Edge::getWeight() { return weight; }

int Edge::getInd1() { return ind1; }

int Edge::getInd2() { return ind2; }

bool cmp(Edge* edge1, Edge* edge2) {
	if (edge1 == nullptr) { return false; } if (edge2 == nullptr) { return true; }
	double weight1 = edge1->getWeight();
	double weight2 = edge2->getWeight();
	return weight1<weight2;
}

bool cmpHeap(indDistPair a, indDistPair b) {
	return a.dist<b.dist;
}

std::string printPair(indDistPair a) {
	std::string s = std::to_string(a.dist);
	return s;
}

vector<int>* createPairing(unsigned int size,unsigned int deg){
    vector<int>* points = new vector<int>(size*deg,-1);
    bool foundPairing=false;
    if(deg==0){
        return points;
    }
    while(!foundPairing){
        srand(time(nullptr));
        bool currPairing=true;
        for(unsigned int i=0;i<size && currPairing;i++){
            for(unsigned int j=0;j<deg && currPairing;j++){
                if((*points)[i*deg+j]==-1){
                    int rndIndex;
                    unsigned int count=0;
                    while(currPairing && !goodPairing(points,i*deg+j,rndIndex = i*deg+rand()%((size-i)*deg),deg)){
                        count++;
                        if(count==(size*deg)/2){ /*Test if a good pairing still exsists*/
                            currPairing = freePairingExists(points,deg,i);
                        }
                    }
                    if(currPairing){
                        (*points)[i*deg+j] = rndIndex;
                        (*points)[rndIndex] = i*deg+j;
                        if(i==size-1 && j==deg-1){
                            foundPairing=true;
                        }
                    }
                }else{
                    if(i==size-1 && j==deg-1){
                        foundPairing=true;
                    }
                }
            }
        }
        if(!foundPairing){
            delete(points);
            points = new vector<int>(size*deg,-1);
        }
    }
    return points;
}

bool goodPairing(vector<int>* pairing,int i,int j,int deg){
    if(j<=i || (*pairing)[j]!=-1){
        return false;
    }
    int ind1 = i/deg;
    int ind2 = j/deg;
    if(ind1==ind2){return false;}
    for(int k=ind1*deg;k<(ind1+1)*deg;k++){
        if((*pairing)[k]!=-1 && ((*pairing)[k])/deg==ind2){
            return false;
        }
    }
    return true;
}

bool freePairingExists(vector<int>* pairing,unsigned int deg,unsigned int ind){
    for(unsigned int i = (ind+1)*deg;i<pairing->size();i++){
        if(goodPairing(pairing,i,(*pairing)[i],deg)){
            return true;
        }
    }
    return false;
}

preEdge translateLine(std::string tString){
    std::stringstream stream;
    stream << tString;
    std::string currSeg;
    int ind1,ind2;
    double weight;
    int i=0;
    while(std::getline(stream,currSeg,',')){
        if(i==0){
            ind1 = std::stoi(currSeg);
        }
        if(i==1){
            ind2 = std::stoi(currSeg);
        }
        if(i==2){
            weight = std::stod(currSeg);
        }
        if(i>2){
            throw std::invalid_argument("Invalid file");
        }
        i++;
    }
    return preEdge(ind1,ind2,weight);
}
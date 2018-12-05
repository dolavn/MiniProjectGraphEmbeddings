#include "Main.h"
#include "Heap.h"
#include "Graph.h"
#include "GraphGenerator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

typedef struct Settings{
    GraphGenerator* generator;
    std::ofstream* out;
    std::ofstream* outSpanner;
    int num;
    std::vector<float> t;
    std::string csv;
    
    Settings(GraphGenerator* generator,std::ofstream* out,std::ofstream* outSpanner,int num,std::vector<float> t,std::string csv):generator(generator),out(out),outSpanner(outSpanner),num(num),t(t),csv(csv){
    }
    
    Settings(const Settings& other):generator(other.generator),out(other.out),outSpanner(other.outSpanner),num(other.num),t(other.t),csv(other.csv){
    }
    
    Settings& operator=(const Settings& other){
        if(this==&other){
            return *this;
        }else{
            clear();
            generator = other.generator;
            out = other.out;
            outSpanner = other.outSpanner;
            num = other.num;
            csv = other.csv;
            t = other.t;
        }
        return *this;
    }
    
    ~Settings(){
        clear();
    }
    
    void clear(){
         if(generator!=nullptr){
            delete(generator);
            generator = nullptr;
        }
        if(out!=nullptr){
            delete(out);
            out = nullptr;
        }
        if(outSpanner!=nullptr){
            delete(outSpanner);
            outSpanner = nullptr;
        }   
    }
} Settings;

enum GraphType{
    RANDOM=1,REGULAR=2,BIPARTITE=3,PLANAR=4,ERROR_TYPE=0
};

bool cmp(int a, int b);
std::string toString(int a);

bool cmp(int a, int b) {
	return a < b;
}

std::string toString(int a) {
	return std::to_string(a);
}

std::vector<float> splitString(std::string str);

void printPointers(std::vector<int*> p) {
	for (unsigned int i = 0; i < p.size(); i++) { std::cout << *p[i] << std::endl; }
}

Settings* getSettings(int argc,char** argv);

int main(int argc, char** argv) {
    Settings* settings;
    try{
        settings = getSettings(argc,argv);
    }catch(std::runtime_error e){
        std::cout << "Error occured:" << e.what() << std::endl;
        return 0;
    }catch(std::invalid_argument e){
        std::cout << "Invalid argument given:" << e.what() << std::endl;
        return 0;
    }
    std::ofstream csvStream;
    std::string csv = settings->csv;
    std::vector<float> tValues = settings->t;
    std::cout << "Generating " << settings->num << " " << settings->generator->toString() << std::endl;
    if(csv!=""){
        csvStream.open(csv);
        csvStream << "Settings:," << settings->generator->toString() << "\n";
        csvStream << "T values:";
        for(unsigned int j=0;j<tValues.size();j++){
            csvStream << "," << settings->t[j];
        }
        
        csvStream << "\nGraph num of edges,Graph weight";
        for(unsigned int j=0;j<tValues.size();j++){
            csvStream << "," << tValues[j] << "-spanner num of edges," << tValues[j] << "-spanner weight";
        }
        csvStream << ",MST num of edges,MST weight\n";
    }
    for(unsigned int j=0;j<tValues.size();j++){
        std::cout << "Creating " << settings->t[j] << "-spanner for each" << std::endl;   
    }
    for(int i=0;i<settings->num;i++){
        Graph* b = settings->generator->generateGraph();
        std::cout << "Graph number " << (i+1) << " generated" << std::endl;
        std::cout << "Graph number of edges:" << b->numOfEdges() << std::endl;
        std::cout << "Graph weight:" << b->getGraphWeight() << std::endl;
        if(settings->out!=nullptr){ 
            b->saveToFile(*settings->out);
        }
        if(csv!=""){
            csvStream << b->numOfEdges() << ",";
            csvStream << b->getGraphWeight();
        }
        for(unsigned int j=0;j<tValues.size();j++){
            Graph* spanner = b->createSpanner(tValues[j]);
            std::cout << "Spanner number of edges:" << spanner->numOfEdges() << std::endl;
            std::cout << "Spanner weight:" << spanner->getGraphWeight() << std::endl;
            if(csv!=""){
                csvStream << "," << spanner->numOfEdges() << "," << spanner->getGraphWeight();
            }
	      if(settings->outSpanner!=nullptr){
                spanner->saveToFile(*settings->outSpanner);
            }
            delete(spanner);
        }
        Graph* mst = b->createMST();
        std::cout << "MST number of edges:" << mst->numOfEdges() << std::endl;
        std::cout << "MST weight:" << mst->getGraphWeight() << std::endl;
        if(csv!=""){
            csvStream << "," << mst->numOfEdges() << "," << mst->getGraphWeight() << "\n";
        }
        delete(b);
        delete(mst);
    }
    if(csv!=""){
        csvStream.close();
    }
    if(settings->out!=nullptr){
        settings->out->close();
    }
    if(settings->outSpanner!=nullptr){
        settings->outSpanner->close();
    }
    return 0;
}

Settings* getSettings(int argc,char** argv){
    GraphType type = ERROR_TYPE;
    int typeInd=-1;
    int num=1; /*Usage -n (num)*/
    std::string csv="";
    std::ofstream* ostream=nullptr;
    std::ofstream* outSpannerStream=nullptr;
    std::string input = "";
    std::vector<float> t = std::vector<float>();
    GraphGenerator* generator=nullptr;
    for(int i=0;i<argc;i++){
        std::string currString = std::string(argv[i]);
        if(currString=="-reg"){
            type = REGULAR;
            typeInd = i;
        }
        if(currString=="-rnd"){
            type = RANDOM;
            typeInd = i;
        }
        if(currString=="-bip"){
            type = BIPARTITE;
            typeInd = i;
        }
        if(currString=="-pln"){
            type = PLANAR;
            typeInd = i;
        }
        if(currString=="-n"){
            num = std::stoi(std::string(argv[i+1]),nullptr,10);
        }
        if(currString=="-t"){
            std::vector<float> currFloat = splitString(std::string(argv[i+1]));
            for(unsigned int j=0;j<currFloat.size();j++){
                t.push_back(currFloat[j]);
            }
        }
        if(currString=="-o"){
            ostream = new std::ofstream(std::string(argv[i+1]));
        }
        if(currString=="-os"){
            outSpannerStream = new std::ofstream(std::string(argv[i+1]));
        }
        if(currString=="-i"){
            input = std::string(argv[i+1]);
        }
        if(currString=="-csv"){
            csv = std::string(argv[i+1]);
        }
    }
    if(num<=0){
        throw std::runtime_error("Illegal number of expirements");
    }
    int nodesNum1,nodesNum2,deg;
    float p,minWeight,maxWeight;
    if(input!=""){
        FileGraphGenerator* fGenerator = new FileGraphGenerator(input);
        num = fGenerator->size();
        generator = fGenerator;
    }
    if(generator==nullptr){
        switch(type){
            case RANDOM: /*Usage: -rnd (nodesNum) (p) (maxWeight)*/
                nodesNum1 = std::stoi(std::string(argv[typeInd+1]));
                p = std::stof(std::string(argv[typeInd+2]));
                minWeight = std::stof(std::string(argv[typeInd+3]));
                maxWeight = std::stof(std::string(argv[typeInd+4]));
                generator = new RandomGraphGenerator(nodesNum1,p,minWeight,maxWeight);
                break;
            case REGULAR: /*Usage: -reg (nodesNum) (deg) (maxWeight)*/
                nodesNum1 = std::stoi(std::string(argv[typeInd+1]));
                deg = std::stoi(std::string(argv[typeInd+2]));
                minWeight = std::stof(std::string(argv[typeInd+3]));
                maxWeight = std::stof(std::string(argv[typeInd+4]));
                generator = new RegularGraphGenerator(nodesNum1,deg,minWeight,maxWeight);
                break;
            case BIPARTITE: /*Usage -bip (nodesNum v1) (nodesNum v2) (p) (maxWeight)*/
                nodesNum1 = std::stoi(std::string(argv[typeInd+1]));
                nodesNum2 = std::stoi(std::string(argv[typeInd+2]));
                p = std::stof(std::string(argv[typeInd+3]));
                minWeight = std::stof(std::string(argv[typeInd+4]));
                maxWeight = std::stof(std::string(argv[typeInd+5]));  
                generator = new BipartiteGraphGenerator(nodesNum1,nodesNum2,p,minWeight,maxWeight);
                break;
            case PLANAR: /*Usage -pln (nodesNum) (p) (minWeight) (maxWeight) */
                nodesNum1 = std::stoi(std::string(argv[typeInd+1]));
                p = std::stof(std::string(argv[typeInd+2]));
                minWeight = std::stof(std::string(argv[typeInd+3]));
                maxWeight = std::stof(std::string(argv[typeInd+4]));
                generator = new PlanarGraphGenerator(nodesNum1,p,minWeight,maxWeight);
                break;
            case ERROR_TYPE: /*No argument given*/
                throw std::runtime_error("No graph type given");
                break;
        }
    }
    return new Settings(generator,ostream,outSpannerStream,num,t,csv);
}

std::vector<float> splitString(std::string tString){
    std::vector<float> ans = std::vector<float>();
    std::stringstream stream;
    stream << tString;
    std::string currSeg;
    while(std::getline(stream,currSeg,',')){
        float currFloat = std::stof(currSeg);
        if(currFloat<1){
            throw std::runtime_error("Invalid t value");
        }
        ans.push_back(std::stof(currSeg));
    }
    return ans;
}

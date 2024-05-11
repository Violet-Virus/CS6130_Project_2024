#ifndef NORMAL_BIPARTITE_GRAPH_H
#define NORMAL_BIPARTITE_GRAPH_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <queue>
#include <limits>
#include <cstring> 
#include "Matching.h"

class NormalBipartiteGraph {
private:
    std::unordered_map<std::string, std::unordered_set<std::string>> adjList;
    std::unordered_set<std::string> setA;
    std::unordered_set<std::string> setB;
    int M = INT16_MAX; // The number of vertices in setA
    int N = INT16_MAX; // The number of vertices in setB

    // Helper function
    int maxBPM(std::vector<std::vector<bool>> bpGraph, std::vector<std::pair<int, int>>& matchedPairs);
    // Helper function
    bool bpm(std::vector<std::vector<bool>> bpGraph, int u, bool seen[], int matchR[], 
        std::vector<std::pair<int, int>>& matchedPairs,
        std::unordered_map<int, int>& jobApplicantMap, std::unordered_map<int, int>& applicantJobMap);
    // Helper function
    void findVertexReachableViaAltPath(std::unordered_set<std::string>& criticalSet, 
        std::unordered_map<std::string, std::string> Matching, std::string vertex);
    // Returns the key given the value
    std::string getKeyByValue(const std::unordered_map<std::string, std::string>& map, const std::string& value);

public:
    // Adds a vertex to the appropriate partition
    void addVertex(const std::string& vertex, bool isSetA); 
    // Adds an edge to the graph
    void addEdge(const std::string& u, const std::string& v); 
    // Removes an edge from the graph
    void removeEdge(const std::string& u, const std::string& v); 

    // Prints the graph
    void printGraph(); 
    // Prints the matching
    void printMatching(const std::unordered_map<std::string, std::string>& matching) const; 

    // Returns the set of neighbours to a vertex
    std::unordered_set<std::string> getNeighbors(const std::string& vertex); 

    // Returns true if the vertex is connected to another vertex
    bool isMatched(const std::string& vertex);
    
    // Computes and returns the MaxMatching
    std::unordered_map<std::string, std::string> computeMaxMatching();
    // Returns boolean value of whether the Graph has a Perfect Matching
    bool hasPerfectMatching();

    // Returns the set of vertices in setA present in the criticalSet
    std::unordered_set<std::string> getCriticalASetVertices(std::unordered_map<std::string, std::string> Matching);
    // Returns the set of vertices in setB which are neighbours of the given set of Vertices
    std::unordered_set<std::string> getNeighbors_OfASet(std::unordered_set<std::string> set);
};


#endif

#include "NormalBipartiteGraph.h"

void NormalBipartiteGraph::addVertex(const std::string& vertex, bool isSetA) 
{
    if (isSetA) 
    {
        setA.insert(vertex);
    } 
    else 
    {
        setB.insert(vertex);
    }
}

void NormalBipartiteGraph::addEdge(const std::string& u, const std::string& v) 
{
    adjList[u].insert(v);
    adjList[v].insert(u);
}

void NormalBipartiteGraph::removeEdge(const std::string& u, const std::string& v) 
{
    adjList[u].erase(v);
    adjList[v].erase(u);
}

void NormalBipartiteGraph::printGraph() 
{
    std::cout << "Vertices in Set A: ";
    for (const std::string& vertex : setA) 
    {
        std::cout << vertex << " ";
    }

    std::cout << "\nVertices in Set B: ";
    for (const std::string& vertex : setB) 
    {
        std::cout << vertex << " ";
    }

    std::cout << "\nEdges:\n";
    for (const auto& [vertex, neighbors] : adjList) 
    {
        std::cout << "Vertex " << vertex << " -> ";
        for (const std::string& neighbor : neighbors) 
        {
            std::cout << neighbor << " ";
        }
        std::cout << std::endl;
    }
}

std::unordered_set<std::string> NormalBipartiteGraph::getNeighbors(const std::string& vertex) 
{
    if (adjList.find(vertex) != adjList.end()) 
    {
        return adjList[vertex];
    }
    return {};
}

bool NormalBipartiteGraph::isMatched(const std::string& vertex) 
{
    if (adjList.find(vertex) != adjList.end() && !adjList[vertex].empty()) 
    {
        return true;
    }
    return false;
}

void NormalBipartiteGraph::printMatching(const std::unordered_map<std::string, std::string>& matching) const 
{    
    std::cout << "Max Matching:\n";
    for (const auto& pair : matching) 
    {
        std::cout << pair.first << " -> " << pair.second << std::endl;
    }
}

bool NormalBipartiteGraph::bpm(std::vector<std::vector<bool>> bpGraph, int u,
         bool seen[], int matchR[], std::vector<std::pair<int, int>>& matchedPairs,
         std::unordered_map<int, int>& setBsetAMap, std::unordered_map<int, int>& setAsetBMap) 
{
    for (int v = 0; v < N; v++) 
    {

        if (bpGraph[u][v] && !seen[v]) 
        {

            seen[v] = true;
            if (matchR[v] < 0 || bpm(bpGraph, matchR[v], seen, matchR, matchedPairs,
                                     setBsetAMap, setAsetBMap)) 
            {
                if (matchR[v] >= 0) 
                {
                    // Remove previous match for setB_element v
                    int prevsetA_element = matchR[v];
                    auto it = matchedPairs.begin();
                    while (it != matchedPairs.end()) 
                    {
                        if (it->second == v && it->first != u)
                        {
                            it = matchedPairs.erase(it);
                        }
                        else
                        {
                            ++it;
                        }
                    }
                    setBsetAMap.erase(v);
                    setAsetBMap.erase(prevsetA_element);
                }

                matchR[v] = u;
                matchedPairs.push_back({u, v}); // Store matched pair
                setBsetAMap[v] = u;
                setAsetBMap[u] = v;
                return true;
            }
        }
    }

    return false;
}


int NormalBipartiteGraph::maxBPM(std::vector<std::vector<bool>> bpGraph, std::vector<std::pair<int, int>>& matchedPairs)
{
    int matchR[N];
    memset(matchR, -1, sizeof(matchR));
    std::unordered_map<int, int> setBsetAMap;    // Maps setB_element to setA_element
    std::unordered_map<int, int> setAsetBMap;    // Maps setA_element to setB_element

    int result = 0;
    for (int u = 0; u < M; u++) 
    {
        bool seen[N];
        memset(seen, 0, sizeof(seen));
        if (bpm(bpGraph, u, seen, matchR, matchedPairs, setBsetAMap, setAsetBMap))
        {
            result++;
        }
    }

    return result;
}

std::unordered_map<std::string, std::string> NormalBipartiteGraph::computeMaxMatching() 
{

    M = setA.size();
    N = setB.size();
    std::vector<std::string> setA_vertices;
    std::vector<std::string> setB_vertices;
    for (const std::string& vertex : setA) 
    {
        setA_vertices.push_back(vertex);
    }
    for (const std::string& vertex : setB) 
    {
        setB_vertices.push_back(vertex);
    }

    std::unordered_map<std::string, std::string> matching;
    std::vector<std::pair<int, int>> matchedPairs;

    std::vector<std::vector<bool>> bpGraph(M, std::vector<bool>(N, false));

    for (int i = 0; i < M; ++i) 
    {
        const std::string& vertexA = setA_vertices[i];
        for (int j = 0; j < N; ++j) {
            const std::string& vertexB = setB_vertices[j];
            bool edgeExists = false;

            auto itA = adjList.find(vertexA);
            if (itA != adjList.end()) 
            {
                const std::unordered_set<std::string>& neighborsA = itA->second;
                edgeExists = neighborsA.find(vertexB) != neighborsA.end();
            } 
            else 
            {
                edgeExists = false;
            }
            bpGraph[i][j] = edgeExists;
        }
    }
    maxBPM(bpGraph, matchedPairs);

    // Convert matchedPairs into matching map
    for (const auto& pair : matchedPairs) 
    {
        matching[setA_vertices[pair.first]] = setB_vertices[pair.second];
    }

    return matching;
}

bool NormalBipartiteGraph::hasPerfectMatching() 
{
    if(setA.size() == setB.size()) 
    {
        if(computeMaxMatching().size() == setA.size()) 
        {
            return true;
        }
    }

    return false;
}

std::string NormalBipartiteGraph::getKeyByValue(const std::unordered_map<std::string, std::string>& map, const std::string& value) 
{
    for (const auto& pair : map) 
    {
        if (pair.second == value) 
        {
            return pair.first; 
        }
    }

    return ""; 
}

void NormalBipartiteGraph::findVertexReachableViaAltPath(std::unordered_set<std::string>& criticalSet, 
    std::unordered_map<std::string, std::string> Matching, std::string vertex) 
{
    std::unordered_set<std::string> vertex_neighbours = getNeighbors(vertex);
    for (const std::string& neigh_vertex : vertex_neighbours) 
    {
        if(neigh_vertex != Matching.find(vertex)->second) 
        {
            auto crtVertex = getKeyByValue(Matching, neigh_vertex);
            criticalSet.insert(crtVertex);
            findVertexReachableViaAltPath(criticalSet, Matching, crtVertex);
        }
    }
}

std::unordered_set<std::string> NormalBipartiteGraph::getCriticalASetVertices(std::unordered_map<std::string, std::string> Matching)
{
    std::unordered_set<std::string> criticalSet;

    for (const std::string& vertex : setA) 
    {
        if (Matching.find(vertex) == Matching.end()) 
        {
            criticalSet.insert(vertex);
        }
    }

    for (const std::string& vertex : criticalSet) 
    {
        std::unordered_set<std::string> vertex_neighbours = getNeighbors(vertex);
        for (const std::string& neigh_vertex : vertex_neighbours) 
        {
            auto crtVertex = getKeyByValue(Matching, neigh_vertex);
            criticalSet.insert(crtVertex);
            findVertexReachableViaAltPath(criticalSet, Matching, crtVertex);
        }
    }

    return criticalSet;
}

std::unordered_set<std::string> NormalBipartiteGraph::getNeighbors_OfASet(std::unordered_set<std::string> set) 
{
    std::unordered_set<std::string> set_neigh;

    for (const std::string& vertex : set) 
    {
        std::unordered_set<std::string> vertex_neighbours = getNeighbors(vertex);
        for (const std::string& neigh_vertex : vertex_neighbours) 
        {
            set_neigh.insert(neigh_vertex);
        }
    }

    return set_neigh;
}






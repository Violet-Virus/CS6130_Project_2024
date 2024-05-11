#include <GraphReader.h>
#include <Utils.h>
#include <Vertex.h>
#include <StronglyStableMatching.h>
#include <iostream>
#include "TestDefs.h"
#include <NormalBipartiteGraph.h>

TEST_CASE("StronglyStableMatching strong_stable", "[matching_stronglystablematching]") {
    auto G = read_graph(get_filepath(get_resources_dir(), "/strong_stable.txt"));
    StronglyStableMatching sm(G);
    auto M = sm.compute_matching();

    if(M.size() != 0)
    {
        std::cout<<sm.verify_strong_stable(M)<<"\n";
    }

}

TEST_CASE("StronglyStableMatching strong_stable_1", "[matching_stronglystablematching]") {
    auto G = read_graph(get_filepath(get_resources_dir(), "/strong_stable_1.txt"));
    StronglyStableMatching sm(G);
    auto M = sm.compute_matching();

    if(M.size() != 0)
    {
        std::cout<<sm.verify_strong_stable(M)<<"\n";
    }
    
}

TEST_CASE("StronglyStableMatching testing_purpose", "[matching_stronglystablematching]") {
    auto G = read_graph(get_filepath(get_resources_dir(), "/testing_purpose.txt"));
    StronglyStableMatching sm(G);
    auto M = sm.compute_matching();

    if(M.size() != 0)
    {
        std::cout<<sm.verify_strong_stable(M)<<"\n";
    }

}
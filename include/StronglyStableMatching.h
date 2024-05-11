#ifndef STRONGLY_STABLE_MATCHING_H
#define STRONGLY_STABLE_MATCHING_H

#include "MatchingAlgorithm.h"
#include "NormalBipartiteGraph.h"

// Strongly stable matching 
// (a, b) is blocking pair in Strongly Stable Matching if
// 1. pref a : b > M(a) and pref b : a >= M(b), or
// 2. pref a : b >= M(a) and pref b : a > M(b)
class StronglyStableMatching : public MatchingAlgorithm {
private:
    // Checks the stop condtion and return a boolean value
    bool check_stop(NormalBipartiteGraph G_sub, std::shared_ptr<BipartiteGraph> G, 
        std::map<VertexPtr, PartnerList> partnerlist_data, bool A_proposing);

public:
    explicit StronglyStableMatching(std::shared_ptr<BipartiteGraph> G, bool A_proposing=true);
    ~StronglyStableMatching() override = default;

    // Given a matching function returns if it is a Strongly Stable Matching or not
    // If not a Strongly Stable Matching then prints the blocking pairs
    bool verify_strong_stable(Matching &M) const;

    // Compute the Strongly Stable Matching if it exists
    // In the case of no Strongly Stable Matching it returns a empty matching
    Matching compute_matching() override;
};

#endif

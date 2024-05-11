#include "StronglyStableMatching.h"
#include "MatchingAlgorithm.h"
#include "VertexBookkeeping.h"
#include "Popular.h"
#include "Utils.h"
#include "Vertex.h"
#include "NormalBipartiteGraph.h"
#include "Matching.h"
#include <iostream>

StronglyStableMatching::StronglyStableMatching(std::shared_ptr<BipartiteGraph> G, bool A_proposing)
    : MatchingAlgorithm(G, A_proposing)
{}

bool StronglyStableMatching::check_stop(NormalBipartiteGraph G_sub, std::shared_ptr<BipartiteGraph> G, std::map<VertexPtr, PartnerList> partnerlist_data, bool is_A_proposing)
{

  const auto& proposing_partition = is_A_proposing ? G->get_A_partition()
                                                     : G->get_B_partition();

  // Checks all the unmatched elements in the proposing_partition if their partner list is empty of not
  // If empty the algo can stop                                             
  for(const auto& element : proposing_partition) 
  {
    if(partnerlist_data[element.second].size() == 0 && G_sub.getNeighbors(element.second.get()->get_id()).size() == 0) 
    {
      return false;
    }
  }
  return true;
}

Matching StronglyStableMatching::compute_matching() 
{

  std::shared_ptr<BipartiteGraph> G = get_graph();
  auto M = Matching(is_A_proposing());

  // List to store unmatched vertices in the proposing partition
  FreeListType free_list;
  std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
  // The possible partners of the vertices in the proposing partition are stored 
  std::map<VertexPtr, PartnerList> partnerlist_data;

  // choose the partitions from which the vertices will propose
  const auto& proposing_partition = is_A_proposing() ? G->get_A_partition()
                                                     : G->get_B_partition();

  const auto& non_proposing_partition = is_A_proposing() ? G->get_B_partition()
                                                         : G->get_A_partition();

  // Initializing the bookkeeping data and free list for vertices in the proposing partition                                                          
  for (auto &it : proposing_partition) 
  {
    auto v = it.second;
    free_list.push(v);
    int pref_list_size = v->get_preference_list().size();
    int residual = non_proposing_partition.size();
    bookkeep_data[v] = VertexBookkeeping(0, pref_list_size, 0, residual);
  }

  // Initializing the partner list data for vertices in the proposing partition
  for (auto &it : proposing_partition) 
  {
    auto v = it.second;
    PartnerList list = PartnerList();
    auto v_prefS_list = v->get_preference_list().get_prefS();
    for (auto &partner_of_v : v_prefS_list) 
    {
      list.add_partner(partner_of_v.vertex, partner_of_v.rank, 0);
    }
    partnerlist_data[v] = list;
  }

  // Initializing the partner list data for vertices in the non proposing partition
  for (auto &it : non_proposing_partition) 
  {
    auto v = it.second;
    PartnerList list = PartnerList();
    auto v_prefS_list = v->get_preference_list().get_prefS();
    for (auto &partner_of_v : v_prefS_list) 
    {
      list.add_partner(partner_of_v.vertex, partner_of_v.rank, 0);
    }
    partnerlist_data[v] = list;
  }

  // Initializing the sub Graph by adding all the vertices
  // Added proposing_partition vertices in the setA
  // Added non_proposing_partition vertices in the setB
  NormalBipartiteGraph G_sub;
  for (auto &it : proposing_partition) 
  {
    auto v = it.first;
    G_sub.addVertex(v, true);
  }
  for (auto &it : non_proposing_partition) 
  {
    auto v = it.first;
    G_sub.addVertex(v, false);
  }

  do 
  {
    // Traverse the list of vertices in the proposing_partition in the setA who are not matched i.e in the free list
    while (not free_list.empty()) 
    {
      // remove the vertex from the free list
      auto u = remove_from_free_list(free_list, bookkeep_data);

      const auto &u_pref_list = u->get_preference_list();
      const auto &u_prefS_list = u_pref_list.get_prefS();

      auto &u_data = bookkeep_data[u];
      auto &partner_list_u = partnerlist_data[u];

      // while the vertex u has not exhausted its preference list and it is not matched
      while (!u_data.is_exhausted() && !G_sub.isMatched(u.get()->get_id())) 
      {

        VertexPtr v;

        // Case where the vertices at the begin of the pref list are tied 
        if (u_pref_list.is_tied(u_data.begin)) 
        {

          auto ties = u_pref_list.get_ties(u_data.begin);
          auto num_ties = ties.size();

          // for all the vertices tied in the list
          for(int idx = 0; idx < num_ties; ++idx) {

            v = ties[idx].vertex;
            // check if vertex is there in the partner list 
            if(partner_list_u.find(v) != partner_list_u.end()) 
            {

              // check if that vertex is matched 

              // if matched check the pref list to see which is the prefered partner
              if(G_sub.isMatched(v.get()->get_id())) 
              {

                auto matchedVertices = G_sub.getNeighbors(v.get()->get_id());
                auto element = *matchedVertices.begin();

                // matchedV is the vertex matched to v
                auto matchedV = G.get()->get_A_partition().find(element);
                auto pref = v.get()->get_preference_list().prefers(u, matchedV->second);

                // if u is better then the matchedV remove all prior engagements to v and get it engaged to u
                if(pref == cBetter)
                {
                  //changes made to the sub graph
                  for (const std::string& element_to_delete : matchedVertices) 
                  {

                    G_sub.removeEdge(v.get()->get_id(), element_to_delete);

                    // if that vertex : element_to_delete is now  not connected to any other vertex then 
                    // it is added back to the free list
                    if(!G_sub.isMatched(element_to_delete))
                    {
                      add_to_free_list(free_list, G.get()->get_A_partition().find(element_to_delete)->second);
                    }
                  }
                  G_sub.addEdge(u.get()->get_id(),v.get()->get_id());

                  //modification of the partner list to delete vertices less prefered than u in v's list
                  auto v_prefS = (v->get_preference_list()).get_prefS();
                  for(auto &partner : v_prefS)
                  {
                    auto partner_vertex = partner.vertex;
                    if(v.get()->get_preference_list().prefers(u, partner_vertex) == cBetter)
                    {
                      auto &partner_list_v = partnerlist_data[v];
                      partner_list_v.remove(partner_vertex);
                      auto &partner_list_partner_vertex = partnerlist_data[partner_vertex];
                      partner_list_partner_vertex.remove(v);
                    }
                  }
                } 
                else if (pref == cEqual) 
                {
                  // if u has the same preference as the matchedV then then edge is added to the sub graph
                  // no modifications made in the partner lists
                  G_sub.addEdge(u.get()->get_id(),v.get()->get_id());
                }
              } 
              else 
              {
                // case where the vertex v is not matched
                // the edge between u and v is added to the graph
                G_sub.addEdge(u.get()->get_id(),v.get()->get_id());
                
                //modification of the partner list to delete vertices less prefered than u in v's list
                auto v_prefS = (v->get_preference_list()).get_prefS();
                for(auto &partner : v_prefS) 
                {
                  auto partner_vertex = partner.vertex;
                  if(v.get()->get_preference_list().prefers(u, partner_vertex) == cBetter) 
                  {
                    auto &partner_list_v = partnerlist_data[v];
                    partner_list_v.remove(partner_vertex);
                    auto &partner_list_partner_vertex = partnerlist_data[partner_vertex];
                    partner_list_partner_vertex.remove(v);
                  }
                }
              }
            }
          }
        } 
        // Case where the vertices at the begin of the pref list is not tied 
        else 
        {
          v = u_pref_list.at(u_data.begin).vertex;
          // check if vertex is there in the partner list  
          if(partner_list_u.find(v) != partner_list_u.end()) 
          {
            // check if that vertex is matched 

            // if matched check the pref list to see which is the prefered partner
            if(G_sub.isMatched(v.get()->get_id())) 
            {
              auto matchedVertices = G_sub.getNeighbors(v.get()->get_id());
              auto element = *matchedVertices.begin();
              
              // matchedV is the vertex matched to v
              auto matchedV = G.get()->get_A_partition().find(element);
              auto pref = v.get()->get_preference_list().prefers(u, matchedV->second);

              // if u is better then the matchedV remove all prior engagements to v and get it engaged to u
              if(pref == cBetter)
              {
                //changes made to the sub graph
                for (const std::string& element_to_delete : matchedVertices) 
                {

                  G_sub.removeEdge(v.get()->get_id(), element_to_delete);

                  // if that vertex : element_to_delete is now  not connected to any other vertex then 
                  // it is added back to the free list
                  if(!G_sub.isMatched(element_to_delete))
                  {
                    // std::cout<<"FOR THE VERTEX "<<element_to_delete<<" is added back to the free list\n";
                    add_to_free_list(free_list, G.get()->get_A_partition().find(element_to_delete)->second);
                  }
                }
                G_sub.addEdge(u.get()->get_id(),v.get()->get_id());

                //modification of the partner list to delete vertices less prefered than u in v's list
                auto v_prefS = (v->get_preference_list()).get_prefS();
                for(auto &partner : v_prefS)
                {
                  auto partner_vertex = partner.vertex;
                  if(v.get()->get_preference_list().prefers(u, partner_vertex) == cBetter)
                  {
                    auto &partner_list_v = partnerlist_data[v];
                    partner_list_v.remove(partner_vertex);
                    auto &partner_list_partner_vertex = partnerlist_data[partner_vertex];
                    partner_list_partner_vertex.remove(v);
                  }
                }
              } 
              else if (pref == cEqual) 
              {
                // if u has the same preference as the matchedV then then edge is added to the sub graph
                // no modifications made in the partner lists
                G_sub.addEdge(u.get()->get_id(),v.get()->get_id());
              }
            } 
            else 
            {
              // case where the vertex v is not matched
              // the edge between u and v is added to the graph
              G_sub.addEdge(u.get()->get_id(),v.get()->get_id());

              //modification of the partner list to delete vertices less prefered than u in v's list
              auto v_prefS = (v->get_preference_list()).get_prefS();
              for(auto &partner : v_prefS)
              {
                auto partner_vertex = partner.vertex;
                if(v.get()->get_preference_list().prefers(u, partner_vertex) == cBetter)
                {
                  auto &partner_list_v = partnerlist_data[v];
                  partner_list_v.remove(partner_vertex);
                  auto &partner_list_partner_vertex = partnerlist_data[partner_vertex];
                  partner_list_partner_vertex.remove(v);
                }
              }
            }
          }
        } 

        //increase the begin index
        u_data.begin += 1;
      }

      // if the vertex u is not matched after exhausting its pref list the graph
      // does not admit an Strongly Stable Matching
      if(!G_sub.isMatched(u.get()->get_id()) && u_data.is_exhausted())
      {
        goto end_of_while;
      }
    }

    // If we get an subgraph will all the vertices on the proposing side have edges 
    // we try to compute a max matching on the sub graph

    // if there is a perfect matching it is returned 
    if(G_sub.hasPerfectMatching())
    {
      // MaxMatching which is a perfect matching is stored in the matching_sub
      std::unordered_map<std::string, std::string> matching_sub = G_sub.computeMaxMatching();

      // Construction of the matching to be returned 
      auto Matching_final = Matching(is_A_proposing());

      for (const auto& pair : matching_sub) 
      {
        auto a_v = G.get()->get_A_partition().find(pair.first)->second;
        auto b_v = G.get()->get_B_partition().find(pair.second)->second;

        auto rank_b_in_a_list = (RankType)a_v.get()->get_preference_list().find_index(b_v);
        auto rank_a_in_b_list = (RankType)b_v.get()->get_preference_list().find_index(a_v);

        Matching_final.add_partner(G.get()->get_A_partition().find(pair.first)->second, G.get()->get_B_partition().find(pair.second)->second, rank_b_in_a_list, 0);
        Matching_final.add_partner(G.get()->get_B_partition().find(pair.second)->second, G.get()->get_A_partition().find(pair.first)->second, rank_a_in_b_list, 0);
      }

      // the matching is returned 
      return Matching_final;
    }
    // else the critical set is computed and the neighbours of the critical set 
    // that are vertices in the non_proposing_partition have their engagements broken off
    // and the vertices that are the least prefered are removed from each vertex in that set
    else
    {
      // criticalZ is the critical set of vertices
      auto criticalZ = G_sub.getCriticalASetVertices(G_sub.computeMaxMatching());

      // criticalZ_neigh are the neighbours of the vertices in the criticalZ set
      auto criticalZ_neigh = G_sub.getNeighbors_OfASet(criticalZ);

      // all engagements to vertices in criticalZ_neigh are broken
      for (const auto& critical_vertex : criticalZ_neigh) 
      {
        auto critical_vertex_neighbours = G_sub.getNeighbors(critical_vertex);
        for (const auto& critical_vertex_neigh : critical_vertex_neighbours) 
        {
          G_sub.removeEdge(critical_vertex, critical_vertex_neigh);
        }
      }

      // all the vertices in the criticalZ_neigh set have their partners 
      // at the tail of the partner list are removed
      for (const auto& critical_vertex : criticalZ_neigh) 
      {
        auto crt_vertex = G.get()->get_B_partition().find(critical_vertex)->second;

        auto &partner_list_critical_vertex = partnerlist_data[crt_vertex];
        auto prefList_critical_vertex = crt_vertex.get()->get_preference_list();
        auto prefS_critical_vertex = prefList_critical_vertex.get_prefS();
        
        auto highest_rank = prefList_critical_vertex.find_index(partner_list_critical_vertex.begin()->vertex);

        // traverse the partner list to find the largest rank i.e the least prefered vertex/vertices rank
        for(const auto& partnerListVertex : partner_list_critical_vertex) 
        {
          auto rank_partnerListVertex = prefList_critical_vertex.find_index(partnerListVertex.vertex);
          if(rank_partnerListVertex > highest_rank) 
          {
            highest_rank = rank_partnerListVertex;
          }
        }

        auto copy_partner_list_critical_vertex = partner_list_critical_vertex;

        // for every vertex pcv in the partner list that is the least prefered is removed from the partner lists
        for(const auto& pcv : copy_partner_list_critical_vertex)
        {
          auto rank_pcv = prefList_critical_vertex.find_index(pcv.vertex);
          if(rank_pcv == highest_rank)
          {
            partner_list_critical_vertex.remove(pcv.vertex);
            partnerlist_data[pcv.vertex].remove(crt_vertex);
          }
        }
      }

      // The free list is updated to add the unmatched vertices in the proposing_partition
      for (auto &it : proposing_partition) 
      {
        auto v = it.second;
        if(G_sub.getNeighbors(v.get()->get_id()).size() == 0)
        {
          free_list.push(v);
        }
      }
    }

  // the while loop stop cond is to check if any vertex in the proposing_partition has exhausted its pref list
  }while(check_stop(G_sub, G, partnerlist_data, is_A_proposing()));  
  
  end_of_while:

  std::cout<<"\n\nTHERE IS NO STRONGLY STABLE MATCHING\n\n";

  return M;
}

//returns true if blocking pair else false
bool check_strong_blocking_pair(VertexPtr a , VertexPtr b , Matching & M) {

  VertexPtr matched_partner_a = nullptr;
  VertexPtr matched_partner_b = nullptr;
  if(M.has_partner(a) == true)
  {
    matched_partner_a = M.get_partner(a);
  }
  if(M.has_partner(b) == true)
  {
    matched_partner_b = M.get_partner(b);
  }

  if(matched_partner_a != nullptr && matched_partner_b != nullptr)
  {
    if(a->get_preference_list().prefers(matched_partner_a,b) == cWorse && b->get_preference_list().prefers(matched_partner_b,a) != cBetter)
    {
      return true;
    }
    else if(a->get_preference_list().prefers(matched_partner_a,b) != cBetter && b->get_preference_list().prefers(matched_partner_b,a) == cWorse)
    {
      return true;
    }
    return false;
  }
  else if(matched_partner_a == nullptr && matched_partner_b != nullptr)
  {
    if(b->get_preference_list().prefers(matched_partner_b,a) != cBetter)
    {
      return true;
    }
    return false;
  }
  else if(matched_partner_a != nullptr && matched_partner_b == nullptr)
  {
    if(a->get_preference_list().prefers(matched_partner_a,b) != cBetter)
    {
      return true;
    }
    return false;
  }
  return true;
}

bool StronglyStableMatching::verify_strong_stable(Matching &M) const
{                                   
  std::shared_ptr<BipartiteGraph> G = get_graph();                                                                                                                           
  const auto& A_partition = G->get_A_partition();                                                    
  std::cout<<"\n\n";
  bool check_StrongSM = true ;
  for (const auto& [_, a] : A_partition) 
  {
    auto a_prefS = (a->get_preference_list()).get_prefS();
    
    for(auto &partner : a_prefS)
    {
      auto b = partner.vertex ;
      
      if(!M.is_matched_to(a,b))
      {
        if(check_strong_blocking_pair(a,b,M))
        {
          check_StrongSM = false;
          std::cout<<"Blocking pair : ("<<a.get()->get_id()<<", "<<b.get()->get_id()<<")"<<std::endl;
        }
      }
    }
  } 
  std::cout<<"\n\n";
  return check_StrongSM ;
}
/****************************************************************************************
 * File: louv.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 26/11/2014

 ***************************************************************************************/

#ifndef SAMOGWAS_LOUVAIN_LOUV_HPP 
#define SAMOGWAS_LOUVAIN_LOUV_HPP

#include <memory>

#include "distance/comparable.hpp"
#include "clustering/clustering.hpp"
#include "clustering/partition.hpp"

#include "community.hpp"
#include "link.hpp"

namespace samogwas
{

namespace louvain {

class MethodLouvain: public AlgoClusteringInterface  {
 public:
  
  typedef SimilarityMatrix Weights; 
  typedef std::shared_ptr<Weights> WeightsPtr;
  typedef std::shared_ptr<Graph> GraphPtr;
  typedef std::shared_ptr<Network> NetworkPtr;
  
  MethodLouvain( WeightsPtr wts);
  virtual Partition run();
  virtual char* name() const;
  virtual void invalidate() {}  
  //protected:
 public:
  /////////
  inline void first_phase();
  /////////
  inline void second_phase();
  /////////
  inline void randomize_order( std::vector<NodeIndex>& vt );
  /////////
  inline void resize_order( std::vector<NodeIndex>& node_order, const size_t sz );
  ////////
  inline CommunityIndex comm_best_gain( 
                                        const NodeIndex& node,
                                        double own_shared,
                                        double& best_gain,
                                        double& best_shared );
  
 protected:  
  // GraphPtr graph;
  NetworkPtr network;
  
  bool changed;
  std::vector<NodeIndex> node_order; // ensures the order of node processing is random
};




inline void MethodLouvain::randomize_order( std::vector<NodeIndex>& vt ) {
  auto engine = std::default_random_engine{};
  std::shuffle(std::begin(vt), std::end(vt), engine);
}



/**
 *
 */
void MethodLouvain::first_phase() {
  //printf("1st_phase starting  --------------- of ------------------- %d\n", (int)network->nbrNodes() );
  printf("\n\n\n---------------------1st_phase starting, modularity: %f---------------------------------------\n--------------------------------------------------------------------------------------------\n", network->modularity());
    
  resize_order(node_order, network->nbrNodes());
  randomize_order(node_order);

  changed = false;
  double local_changed = true;  
  while (local_changed) {
    // printf("1st-1\n");
    local_changed = false;
    for ( auto node: node_order ) {
      // printf("1st-%d-2\n", node);
      double best_gain, best_shared;
      CommunityIndex own_comm = network->communityOf(node);
      // printf("1st-%d-3\n", node);
      double own_shared = network->sharedWeights(node, network->communityOf(node));
      // printf("1st-%d-4\n", node);
      CommunityIndex best_comm = comm_best_gain( node, own_shared, best_gain, best_shared);
      // printf("1st-%d-5\n", node);
      if ( own_comm != best_comm )
      {
        double tmp = network->modularity();
        changed = true;
        local_changed = true;
        network->moveNode( node, best_comm, own_shared, best_shared );
        // printf("moving %d --> %d, modularity: %f -> %f gain: %f\n", node, best_comm, tmp, network->modularity(), best_gain);
      }      
    }    
  }
}


/**
 *
 */
void MethodLouvain::second_phase() {

  // std::cout << "nbrComm: " << network->nbrCommunities() << std::endl;
  auto links = std::make_shared<WeightedLink>(network->nbrCommunities());
  std::vector<CommunityIndex> community_map;
  // for ( auto comm: network->communities() ) printf("quytquyt of %d - count: %d\n", comm, network->community_member_counts[comm]);
  std::copy( network->communities().begin(), network->communities().end(), std::back_inserter(community_map));

  for ( CommunityIndex c_i = 0; c_i < network->nbrCommunities(); ++c_i ) {
    for ( CommunityIndex c_j = 0; c_j < network->nbrCommunities(); ++c_j  ) {
      CommunityIndex old_c_i = community_map[c_i], old_c_j = community_map[c_j];
      links->setWeight(c_i, c_j, network->interCommunityWeight(old_c_i, old_c_j));
    }
  }
  
  auto g = std::make_shared<Graph>( links );
  network->setGraph( g );

  for ( auto vp = g->allLinks(); vp.first != vp.second; ++vp.first ) {
    auto l = g->source(*vp.first), r = g->target(*vp.first);

  }
  printf("2nd_phase done: modularity = %f\n\n\n\n", network->modularity());

}

inline void MethodLouvain::resize_order( std::vector<NodeIndex>& node_order, const size_t sz ) {  
  node_order.resize(sz);
  for ( int i = 0; i < sz; ++ i ) {
    node_order[i] = i;
  }
}

////////////
inline CommunityIndex MethodLouvain::comm_best_gain( 
                                                     const NodeIndex& node,
                                                     double own_shared,
                                                     double& best_gain,
                                                     double& best_shared ) {
  // @todo:  # only consider neighbors of different communities
  double loss = network->modularityLoss( node, own_shared);
  best_shared = -2.0; best_gain = -loss;
  CommunityIndex best_comm = network->communityOf(node), own_comm = network->communityOf(node);

  std::set<CommunityIndex> processed_neigh_comms { best_comm };  
  for ( auto vp = network->adjacentNodes(node); vp.first != vp.second; ++vp.first ) {
    NodeIndex neigh_node = *vp.first;
    CommunityIndex neigh_comm = network->communityOf(neigh_node);
    if ( processed_neigh_comms.find(neigh_comm) != processed_neigh_comms.end() )  continue; // continues if already processed
    processed_neigh_comms.insert(neigh_comm); // already seen it
    
    double shared_weights = network->sharedWeights(node, neigh_comm);
    double gain = network->modularityGain(node, neigh_comm, shared_weights);

    if ( gain > best_gain ) {
      best_gain = gain;
      best_shared = shared_weights;
      best_comm = neigh_comm;        
    }
  }

  return best_comm;
}



}

} // namespace samogawsends here. samogaws


/********************************************************/
#endif // SAMOGAWS_LOUV_HPP







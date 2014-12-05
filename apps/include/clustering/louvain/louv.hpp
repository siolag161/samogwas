/****************************************************************************************
 * File: louv.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 26/11/2014

 ***************************************************************************************/

#ifndef SAMOGWAS_LOUVAIN_LOUV_HPP 
#define SAMOGWAS_LOUVAIN_LOUV_HPP

#include <memory>
#include <map>

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
  typedef std::shared_ptr<Partition> PartitionPtr;
  
  MethodLouvain( WeightsPtr wts);
  virtual Partition run();
  virtual char* name() const;
  virtual void invalidate() {}  
  //protected:
 public:
  /////////
  inline bool first_phase();
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
  
  inline void update_partition( Partition& partition,
                                const Network& ntw,
                                const std::map<CommunityIndex,CommunityIndex>& comm_map );
  
  GraphPtr graph;
  NetworkPtr network;
  
  PartitionPtr partition;
  // std::map<CommunityIndex, CommunityIndex> local_global_map;  
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
bool MethodLouvain::first_phase() {
  
  printf("----------------haha 1s phase - we have %d over %d: %f---------------\n",
         network->modularity(), network->nbrNodes(), network->nbrCommunities() );
  resize_order(node_order, network->nbrNodes());
  changed = false;
  double local_changed = true;  
  while (local_changed) {
    local_changed = false;
    for ( auto node: node_order ) {
      double best_gain, best_shared;
      CommunityIndex own_comm = network->getCommunity(node);
      double own_shared = network->sharedWeights(node, network->getCommunity(node));
      CommunityIndex best_comm = comm_best_gain( node, own_shared, best_gain, best_shared);
      if ( own_comm != best_comm )
      {
        double tmp = network->modularity();
        changed = true;
        local_changed = true;
        network->moveNode( node, best_comm, own_shared, best_shared );
      } 
    }    
  }
  printf("\n----------------hehe end 1st phase -  now %d vars -> %d communitities, with %f modularity ---------------\n\n", network->nbrNodes(), network->nbrCommunities(), network->modularity() );
  return changed;
}


/**
 *
 */
void MethodLouvain::second_phase() {
  auto nbrComms = network->nbrCommunities();
  std::vector<CommunityIndex> new2old_comm_map;
  std::copy( network->communities().begin(), network->communities().end(), std::back_inserter(new2old_comm_map) );
  std::map<CommunityIndex, CommunityIndex> old2new_comm_map; 
 

  for ( CommunityIndex new_comm = 0; new_comm < new2old_comm_map.size(); ++new_comm ) {
    auto old_comm = new2old_comm_map[new_comm];
    old2new_comm_map[old_comm] = new_comm;
  }

  auto links = std::make_shared<WeightedLink>(nbrComms);

  // #pragma omp parallel for 
 
  
  for ( auto vp = network->allLinks(); vp.first != vp.second; ++vp.first ) {
    NodeIndex l = network->source(*vp.first), r = network->target(*vp.first);
    CommunityIndex c_l = old2new_comm_map[network->getCommunity(l)], c_r = old2new_comm_map[network->getCommunity(r)];
    if ( c_l != c_r ) {
      double tmp= links->compute(c_l, c_r) ;
      double curr_w = links->compute(c_l, c_r) + network->weight(l,r);
      links->setWeight(c_l, c_r, curr_w);
    }     
    // 
  }

  for ( CommunityIndex c_i = 0; c_i < nbrComms; ++c_i ) {
    auto old_comm = new2old_comm_map[c_i];
    links->setWeight(c_i, c_i, network->innerCommunityWeight(old_comm) );
    // printf("self_node[%d]: %f\n", c_i, links->compute(c_i,c_i));
    // for ( CommunityIndex c_j = c_i+1; c_j < nbrComms; ++c_j ) {      
    //   printf("weight[%d-%d]: %f\n", c_i, c_j, links->compute(c_i,c_j));
    // }
  }  
  

  update_partition( *partition, *network, old2new_comm_map); 
  auto g = std::make_shared<Graph>(links);
  network->setGraph(g);

  // printf("\n----------------hehe ends 2nd phase----we have %d over %d with %f modularity---------------\n",
  //        network->nbrNodes(), network->nbrCommunities(), network->modularity() );

}

inline void MethodLouvain::resize_order( std::vector<NodeIndex>& node_order, const size_t sz ) {  
  node_order.resize(sz, -1);
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
  CommunityIndex best_comm = network->getCommunity(node), own_comm = network->getCommunity(node);

  std::set<CommunityIndex> processed_neigh_comms { best_comm };  
  for ( auto vp = network->adjacentNodes(node); vp.first != vp.second; ++vp.first ) {
    NodeIndex neigh_node = *vp.first;
    CommunityIndex neigh_comm = network->getCommunity(neigh_node);
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


inline void MethodLouvain::update_partition( Partition& partition,
                                             const Network& ntw,
                                             const std::map<CommunityIndex,CommunityIndex>& old2new_comm_map ) {
  for ( NodeIndex node = 0; node < partition.nbrItems(); ++node ) {
    auto lab = partition.getLabel(node);
    auto old_comm = ntw.getCommunity(lab);
    auto new_comm = old2new_comm_map.at(old_comm);
    partition.setLabel(node, new_comm);    
  }
}


}} // namespace samogawsends here. samogaws


/********************************************************/
#endif // SAMOGAWS_LOUV_HPP









#include <memory>
 
#include "clustering/clustering.hpp" // AlgoClust
#include "clustering/partition.hpp"  // Partition, Clustering
#include "distance/comparable.hpp"

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"

#include "clustering/louvain/graph.hpp"

namespace samogwas
{
namespace louvain {
  
Graph::Graph( WeightMatPtr l, bool keep_self_loops): weights(l), total_weights(0.0), nbr_links(0) {
  initialize(keep_self_loops);
  
}

void Graph::initialize(bool keep_self_loops) {
  
  linked_weights = std::make_shared<Weights>(nbrNodes(), -1);
  self_loops = std::make_shared<Weights>(nbrNodes(), 0.0);
  
  for ( int i = 0; i < this->nbrNodes(); ++i) {      
    addNode();      
  }
  
  for ( int i = 0; i < this->nbrNodes(); ++i) {
    if (keep_self_loops) {
      (*self_loops)[i] = this->weights->compute(i,i);
    }
    total_weights += selfLoopWeight(i)/2;

    for ( int j = i+1; j < this->nbrNodes(); ++j ) {
      Weight w = this->weights->compute(i,j);
      if (w > 0) { 
        addLink(i,j,w);
        total_weights += w;
        nbr_links++;
      }
    }
  }

}
  

/*
 */
Weight Graph::weight( const NodeIndex& i, const NodeIndex& j) const {
  return weights->compute(i,j);
  //return linkWeightMap[boost::edge(i,j,*this).first];
}

/**
 *
 */
Graph::AdjIteRng Graph::adjacentNodes( const NodeIndex& i) const {
  return boost::adjacent_vertices( i, *this );
}


Graph::LinkIteRng Graph::allLinks() const {
  return boost::edges(*this);
}

Graph::OutLinkIteRng Graph::linksFrom( const NodeIndex& i ) const {
  return boost::out_edges(i, *this);
}


/**
 *
 */
size_t Graph::nbrLinks() const {
  return nbr_links;
}

Weight Graph::totalWeights() {
  if (total_weights <= 0) {
    total_weights = 0.0;
    for ( auto vp = allLinks(); vp.first != vp.second; ++vp.first ) {
      total_weights += 2*weight(*vp.first); //linkWeightMap[*vp.first];
    }

    for ( int i = 0; i < this->nbrNodes(); ++i)
      total_weights += selfLoopWeight(i) / 2;
  }  
  return total_weights;
}



}} // namespace

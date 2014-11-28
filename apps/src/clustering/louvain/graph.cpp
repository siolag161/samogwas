
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
  
Graph::Graph( WeightsPtr l): weights(l), total_weights(0.0) {
  for ( int i = 0; i < this->nbrNodes(); ++i)
    boost::add_vertex(*this); // add nodes

  for ( int i = 0; i < this->nbrNodes(); ++i) {
    total_weights += selfLoopWeight(i) / 2;

    for ( int j = i+1; j < this->nbrNodes(); ++j ) {
      double w = this->weights->compute(i,j);
      if (w > 0) { 
        boost::add_edge( i, j, w, *this );
        total_weights += w;
      }
    }

  }

  // linkWeightMap = get(boost::edge_weight, *this);
  nbr_links = boost::num_edges(*this);
  linked_weights.resize(nbrNodes(), -1);
}

/*
 */
double Graph::weight( const NodeIndex& i, const NodeIndex& j) const {
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

double Graph::totalWeights() {
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

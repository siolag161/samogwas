/****************************************************************************************
 * File: graph.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 19/11/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_LOUVAIN_GRAPH_HPP 
#define SAMOGWAS_LOUVAIN_GRAPH_HPP

#include <memory>

// #include "clustering/clustering.hpp" // AlgoClust
// #include "clustering/partition.hpp"  // Partition, Clustering
#include "distance/comparable.hpp"

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"

namespace samogwas
{

namespace louvain {

typedef int NodeIndex;

typedef boost::property<boost::edge_weight_t, float> LinkWeightProperty;
typedef std::pair<NodeIndex, NodeIndex> Link;

typedef boost::adjacency_list<  // adjacency_list is a template depending on :
  boost::listS,               //  The container used for egdes : here, std::list.
  boost::vecS,                //  The container used for vertices: here, std::vector.
  boost::undirectedS,           //  directed or undirected edges ?.
  boost::no_property,                     //  The type that describes a Node.
  LinkWeightProperty               //  The type that describes an Link
  > BaseGraph;

typedef BaseGraph::edge_descriptor LinkIndex;


/**
 *
 */
class Graph: public BaseGraph {

  static const int INVALID_WEIGHT = -1;
  typedef SimilarityMatrix Weights; 
  typedef std::shared_ptr<Weights> WeightsPtr;

  typedef boost::property_map<Graph, boost::edge_weight_t>::type LinkWeightMap;
  
 public:
  Graph( WeightsPtr l);

  /**
   *
   */
  size_t nbrNodes() const { return weights->nbrVariables(); }


  /**
   *
   */
  size_t nbrLinks() const ;

  /**
   *
   */
  double weight( const NodeIndex& i, const NodeIndex& j) const;

  /**
   *
   */
  double weight( const LinkIndex& j) const { return weight( source(j), target(j) ); }


  /**
   *
   */
  double selfLoopWeight( const NodeIndex& i ) const { return weight(i,i); }

    
  double totalWeights();
  
  typedef adjacency_iterator AdjIte;
  typedef std::pair<AdjIte,AdjIte> AdjIteRng; // iterator range for accessing the vertices
  AdjIteRng adjacentNodes( const NodeIndex& i ) const;

  
  NodeIndex source(LinkIndex link) const { return boost::source(link, *this); }
  NodeIndex target(LinkIndex link) const { return boost::target(link, *this); }
    
  typedef edge_iterator LinkIte;
  typedef std::pair<LinkIte, LinkIte> LinkIteRng; // iterator range for accessing the vertices
  LinkIteRng allLinks() const;
  
  typedef out_edge_iterator OutLinkIte;
  typedef std::pair<OutLinkIte, OutLinkIte> OutLinkIteRng; // iterator range for accessing the vertices
  OutLinkIteRng linksFrom( const NodeIndex& i ) const;

  double linkedWeights( const NodeIndex& node) {
    if ( linked_weights[node] == INVALID_WEIGHT ) {
      linked_weights[node] = 0;
      for ( auto vp = linksFrom(node); vp.first != vp.second; ++vp.first ) {
        linked_weights[node] += weight(*vp.first);
      }
      linked_weights[node] += selfLoopWeight(node); // self-loops 
    }
    return linked_weights[node];
  }


  // protected:
  // double computeTotalWeights();
  
 protected:
  /**
   *
   */
  WeightsPtr weights; // gives the weighted link between vertices links[a][b]
  size_t nbr_links;  
  // LinkWeightMap linkWeightMap;
  double total_weights;

  std::vector<double> linked_weights; // weighted_degrees
};


typedef boost::graph_traits<Graph>::vertex_descriptor node_desc;
typedef boost::graph_traits<Graph>::edge_descriptor link_desc;


} // namespace louvain ends here

} // namespace samogwas ends here. samogwas

/****************************************************************************************/


/****************************************************************************************/
#endif // SAMOGWAS_LOUVAIN_GRAPH_HPP

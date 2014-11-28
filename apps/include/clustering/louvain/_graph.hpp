// /****************************************************************************************
//  * File: graph.hpp
//  * Description: 
//  * @author: siolag161 (thanh.phan@outlook.com)
//  * @date: 19/11/2014

//  ***************************************************************************************/
// #ifndef SAMOGWAS_LOUVAIN_GRAPH_HPP 
// #define SAMOGWAS_LOUVAIN_GRAPH_HPP

// #include <memory>

// #include "clustering/clustering.hpp" // AlgoClust
// #include "clustering/partition.hpp"  // Partition, Clustering
// #include "distance/comparable.hpp"
// namespace samogwas
// {

// namespace louvain {

// //template<class Links>
// typedef Index NodeIndex;

// class Graph {
  
//   typedef SimilarityMatrix Weights;
//   typedef Weights* WeightsPtr;
  
//  public:
//   Graph( WeightsPtr l): weights(l) {}

//   /**
//    *
//   */
//   size_t nbrNodes() const { return weights->nbrVariables(); }

//   /**
//    *
//    */
//   //size_t nbrVertices();
  
//   /**
//    *
//    */
//   double weight( const NodeIndex& i, const NodeIndex& j) const;

//   /**
//    *
//    */
//   std::vector<NodeIndex> adjacentNodes( const NodeIndex& i) const;

//  protected:
//   /**
//    *
//    */
//   WeightsPtr weights; // gives the weighted link between vertices links[a][b]
  
// };

// } // namespace louvain ends here

// } // namespace samogwas ends here. samogwas

// /****************************************************************************************/


// /****************************************************************************************/
// #endif // SAMOGWAS_LOUVAIN_GRAPH_HPP

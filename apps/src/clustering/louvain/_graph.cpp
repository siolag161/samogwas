// /****************************************************************************************
//  * File: graph.cpp
//  * Description: 
//  * @author: siolag161 (thanh.phan@outlook.com)   
//  * @date: 19/11/2014

//  ***************************************************************************************/
// #ifndef SAMOGWAS_GRAPH_HPP
// #define SAMOGWAS_GRAPH_HPP

// #include <memory>
 
// #include "clustering/clustering.hpp" // AlgoClust
// #include "clustering/partition.hpp"  // Partition, Clustering
// #include "distance/comparable.hpp"

// #include "clustering/louvain/graph.hpp"

// namespace samogwas 
// {

// namespace louvain {


// /**
//  */
// double Graph::weight( const NodeIndex& i, const NodeIndex& j) const {
//   // std::cout << "size: "  << weights->nbrVariables() << std::endl;
//   // printf("val: %d,%d :%f\n", i, j, weights->compute(i,j));
//   return weights->compute(i,j);
// }

// /**
//  *
//  */
// std::vector<NodeIndex> Graph::adjacentNodes( const NodeIndex& node ) const {
//   std::vector<NodeIndex> adj_nodes;  
//   for ( int i = 0; i < weights->nbrVariables(); ++i) {
//     if ( i != node && (weights->compute(i,node) > 0) ) {
//       adj_nodes.push_back(i);
//     }
//   }
//   return adj_nodes;
// }

// } // namespace louvain ends here

// } // namespace samogwas ends here. samogwas


// /****************************************************************************************/
// #endif // SAMOGWAS_GRAPH_HPP









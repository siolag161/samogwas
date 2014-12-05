/****************************************************************************************
 * File: louv.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 26/11/2014

 ***************************************************************************************/

#include <memory>
#include <algorithm>    // std::random_shuffle
#include <random>

#include "clustering/clustering.hpp"

#include "clustering/louvain/graph.hpp"
#include "clustering/louvain/community.hpp"
#include "clustering/louvain/louv.hpp"

namespace samogwas
{

namespace louvain {

MethodLouvain::MethodLouvain(WeightsPtr wt): changed(true) {
  graph = std::make_shared<Graph>(wt);
  network = std::make_shared<Network>(graph);
  partition = std::make_shared<Partition>();
  for ( NodeIndex n = 0; n < network->nbrNodes(); ++n ) {
    partition->setLabel(n, network->getCommunity(n));
  }
    
}

Partition MethodLouvain::run() {
  // Network partition(network->nbrNodes());
  // Network network(graph);
  while (changed) {
    first_phase();
    if (changed)
      second_phase();
  }

   printf("\n-------------- DONE clustering -------------\n\n");

  // Partition;
  return *partition;
}

char* MethodLouvain::name() const {
  char* name = new char[80];
  sprintf( name, "LOUVAIN%s", "");
  return name;
}


}




} // namespace samogawsends here. samogaws

/****************************************************************************************/

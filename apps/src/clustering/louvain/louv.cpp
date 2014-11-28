
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
  auto g = std::make_shared<Graph>(wt);
  network = std::make_shared<Network>(g);
}

Partition MethodLouvain::run() {
  // Network network(graph);
  while (changed) {
    first_phase();
    second_phase();
  }

  return *network;
}

char* MethodLouvain::name() const {
  char* name = new char[80];
  sprintf( name, "louvain%s", "");
  return name;
}


}




} // namespace samogawsends here. samogaws

/****************************************************************************************/

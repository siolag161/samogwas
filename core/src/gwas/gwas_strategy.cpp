/****************************************************************************************
 * File: gwas_basic_strategy.cpp
 * Description: @todo
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 03/04/2014
 ****************************************************************************************/

#include <gwas/gwas_basic_strategy.hpp>
#include <fltm/graph.hpp>
#include <map>

namespace samogwas {

GWAS_Strategy::Level2Vertices GWAS_Strategy::levels2Vertices( const Graph& graph) {
  Level2Vertices result; result.reserve(10); // 10 is safe enough, hardly above this number of levels
  vertex_iterator vi, vi_end;
  for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {    
    int vertex = *vi;
    result.resize(graph[vertex].level + 1); // modifies the size to be able to accommodate
    result[graph[vertex].level].push_back(vertex);
  }  

  return result;
}


/////////////////////////////////////////////////////////////////////////////////
GWAS_Strategy::Latent2Children GWAS_Strategy::latent2Children(const Graph& graph) {
  std::map<int, int> edgeMap;
  std::map< int, std::vector<int> > latent2children;
  edge_iterator ei, ei_end;

  for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei) {
    edgeMap[ boost::target(*ei,graph) ] = boost::source(*ei, graph);
    latent2children[boost::source(*ei, graph)].push_back(boost::target(*ei, graph));
  }

  return latent2children;
}

}

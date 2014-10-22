/***************************************************************************************
 * File: gwas_basic_strategy.cpp
 * Description: @todo
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 03/04/2014
 ****************************************************************************************/

#include <gwas/gwas_basic_strategy.hpp>
#include <fltm/graph.hpp>
#include <map>

namespace samogwas {

void GWAS_Basic_Strategy::execute( FLTM_Result &result, Matrix& genotype,
                                   Vector &phenotype, stats::StatTest *statTest ) {
  
  const Graph& graph = result.graph;
  auto latent2children = latent2Children(graph);
  edge_iterator ei, ei_end; vertex_iterator vi, vi_end;
  std::map<vertex_t,vertex_t> edgeMap;
  for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei) {
    edgeMap[ boost::target(*ei,graph) ] = boost::source(*ei, graph);
    latent2children[boost::source(*ei, graph)].push_back(boost::target(*ei, graph));
  }

  
  
}


}


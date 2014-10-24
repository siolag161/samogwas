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

void GWAS_Basic_Strategy::execute( Graph& graph,
                                   Matrix& genoMatrix,
                                   Vector& phenoVec ) {
  
  auto latent2children = latent2Children(graph);
  edge_iterator ei, ei_end; vertex_iterator vi, vi_end;
  std::map<vertex_t,vertex_t> edgeMap;
  for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei) {
    edgeMap[ boost::target(*ei,graph) ] = boost::source(*ei, graph);
    latent2children[boost::source(*ei, graph)].push_back(boost::target(*ei, graph));
  }

  // boost::queue<vertex_t> q; 
  // typedef graph_traits<Graph>::vertex_descriptor vertex_t;
  // typedef std::map<vertex_t, boost::default_color_type> ColorMap;
  // typedef boost::associative_property_map<ColorMap> Color;
  // Graph g = getGraph();
  // ColorMap cmap; Color colorMap(cmap);

  
  // std::vector<double> levelThres{ 0.05, 0.1, 0.1 };
  // std::vector<double> scores{ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.01, 0.01};
  // std::shared_ptr<LevelScoreNodeCriterion> criterion_1(new LevelScoreNodeCriterion(levelThres, scores));


  // GWAS_Basic_Visitor( std::shared_ptr<Criterion> eval,
  //                     std::shared_ptr<ScoreMap> scoreMap,
  //                     std::shared_ptr<Color> col)

        
  // boost::breadth_first_visit( g, 10, boost::queue<vertex_t>(), vis, colorMap);

  
}


}


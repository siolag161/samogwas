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

  // std::vector<double> levelThres{ 0.05, 0.1, 0.1 };
  // std::vector<double> scores{ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.01, 0.01};
  // std::shared_ptr<LevelScoreNodeCriterion> criterion_1(new LevelScoreNodeCriterion(levelThres, scores));

  // bfs_visitor vis(criterion_1, colorMap);//(&reached);
  int nbrLevels = latent2children.size();
  boost::queue<vertex_t> q; 

  // for ( auto& vertex: latent2children[nbrLevels-1] ) {    
  //   boost::breadth_first_visit( graph, vertex, q, *visitor, *visitor->color);
  // }

  // GWAS_Basic_Visitor( std::shared_ptr<Criterion> eval, Color& col): evaluator(eval), color(col) {}  

  // const int phenoCard = 2;
  // GWASAssociationTest associationTest( genoMatrix, phenoVec, phenoCard );

  
  // StatTestNodeCriterion statTest( std::vector<double>& thresholdsByLevel,
  //                        std::shared_ptr<stats::GWASAssociationTest> test )


  
}


}


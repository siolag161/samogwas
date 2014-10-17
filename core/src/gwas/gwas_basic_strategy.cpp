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
  vertex_iterator vi, vi_end;
  // size_t max_level = result.level2LatentVars.size();
  auto latent2children = latent2Children(graph);
  // edge_iterator ei, ei_end;

  // for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei) {
  //   edgeMap[ boost::target(*ei,graph) ] = boost::source(*ei, graph);
  //   latent2children[boost::source(*ei, graph)].push_back(boost::target(*ei, graph));
  // }
  //         std::cout << std::endl;

  std::cout << "max_level: " << max_level << std::endl;
  std::vector<int> count(10,0);
  /*        int cand_level = 1;
            std::vector<int> candidates;
            for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
            int vertex = *vi;
            if ( graph[vertex].level == cand_level ) candidates.push_back(vertex);
            if ( edgeMap.find(vertex) == edgeMap.end() ) edgeMap[vertex] = -1;
            count[graph[vertex].level]++;
            }*/

  //         for ( int i = 0; i < 10; ++i) if (count[i]) std::cout << count[i] << std::endl;

  //         std::cout << "-----------------------------------------------------------------------------" << std::endl;
  //         std::vector<int> next_candidates = candidates;
  //         for ( int level = cand_level; level >= 0; --level ) {
  //             std::cout << "current_level: " << level << ". candidates: " << candidates.size() << std::endl;

  //             std::vector<double> dist;
  //             std::vector<double> result;
  //             if ( candidates.size() )
  //                 performTesting( dist, result, statTest, candidates, genoMat, phenotype, graph, perm );

  //             next_candidates.clear();
  //             for( int i = 0; i < candidates.size(); ++i ) {
  //                 int id = candidates[i];
  //                 if ( result[2*i+1] <= thresholds[level] ) {
  //                     os << level << "," << id << ","
  //                             << graph[id].label << ","
  //                             << edgeMap[id] << "," << graph[id].position << ","
  //                             << result[2*i] << "," << result[2*i+1] << std::endl;
  //                     for ( auto& child: latent2children[id] )  {
  //                         next_candidates.push_back(child);
  //                     }
  //                 }
  //             }
  //             candidates = next_candidates;
  //             std::cout << "-----------------------------------------------------------------------------" << std::endl;
  //         }

  //         for ( int level = cand_level + 1; level <= max_level; ++level ) {
  //             candidates.clear();
  //             for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
  //                 int vertex = *vi;
  //                 if ( graph[vertex].level == level ) candidates.push_back(vertex);
  //             }

  //             std::cout << "current_level: " << level << ". candidates: " << candidates.size() << std::endl;
  //             std::vector<double> dist;
  //             std::vector<double> result;
  //             if ( candidates.size() )
  //                 performTesting( dist, result, statTest, candidates, genoMat, phenotype, graph, perm );

  //             next_candidates.clear();
  //             for( int i = 0; i < candidates.size(); ++i ) {
  //                 int id = candidates[i];
  //                 // if ( result[2*i+1] <= thresholds[level] ) {
  //                 os << level << "," << id << ","
  //                         << graph[id].label << ","
  //                         << edgeMap[id] << "," << graph[id].position << ","
  //                         << result[2*i] << "," << result[2*i+1] << std::endl;

  //                 std::cout << level << "," << id << ","
  //                         << graph[id].label << ","
  //                         << edgeMap[id] << "," << graph[id].position << ","
  //                         << result[2*i] << "," << result[2*i+1] << std::endl;

  //                 // }
  //             }
  //         }
  //         os.close();
};
}


// /****************************************************************************************
//  * File: fltm.cpp
//  * Description: 
//  * @author: siolag161 (thanh.phan@outlook.com)
//  * @date: 09/07/2014
//  ***************************************************************************************/

#include <iostream>
#include <thread>
#include <chrono>
#include <memory>
#include <limits>
#include "utils/csv_parser.hpp"

#include <fstream>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> // to obtain the program's name

#include "distance/comparable.hpp"

#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "clustering/louvain/louv.hpp"

#include "gwas/gwas_basic_strategy.hpp"

#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"
#include "fltm/fltm.hpp"

#include "utils/option_printer.hpp"
#include "utils/custom_option_desc.hpp"

#include "utils/logs_utils.hpp"
#include "utils/timer_utils.hpp"

#include "data_load.hpp"
#include "app_options.hpp"

#include "statistics/permutation_test.hpp"
#include "gwas/gwas_basic_strategy.hpp"
#include "gwas/score_node_criteria.hpp"

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <memory>
#include <boost/lockfree/queue.hpp>

#include "gwas_common.hpp"
#include "statistics/test_statistics.hpp"

using namespace utility;
using namespace samogwas;
using namespace stats;

typedef std::vector<size_t> Candidates;
typedef std::vector<Candidates> CandidatesByLevel;
typedef std::vector<int> Cardinalities;

CandidatesByLevel retrieveCandidates( const Graph& g );
Cardinalities retrieveCardinalities( const Graph& g );
Cardinalities retrieveCardinalities( const Graph& g, const Candidates& candidates );


void performTest( const FLTM_Data& fltmData, const Graph& g, std::map<std::string, std::string>& snp2rs,
                  const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                  const CandidatesByLevel& candidates,
                  const Cardinalities &cardinalities,
                  const int permu,
                  const int chr,
                  const std::vector<double> thresholds,
                  boost::filesystem::path outFile );

void performFilteredPermutationTest( const FLTM_Data& fltmData, const Graph& graph,
                                     const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                                     const CandidatesByLevel& candidatesByLevel,
                                     const Cardinalities &cardinalities,
                                     const int permutations,
                                     const int chr,
                                     const std::vector<double> thresholds,
                                     boost::filesystem::path outDir);


void performHasGoodParentPermutationTest( const FLTM_Data& fltmData, const Graph& graph,
                                          const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                                          const CandidatesByLevel& candidatesByLevel,
                                          const Cardinalities &cardinalities,
                                          const int permutations,
                                          const int chr,
                                          const std::vector<double> thresholds,
                                          boost::filesystem::path outDir);

void snpMapping(std::map<std::string, std::string>& snp2rs, std::string& mapFile);

//////////////////////////////////////////////////////////////
int main( int argc, char** argv ) {

  auto pos = getGwasProgramOptions( argc, argv );
  std::cout << "Loading graph data...\n" << std::endl;
  Graph graph;
  BayesGraphLoad graphLoad;
  graphLoad( graph,
             pos.inputLabelFile,
             pos.bayesVertices,
             pos.bayesDist );

  printf("Done loading graph of %d edges and %d vertices\n", boost::num_edges(graph),
         boost::num_vertices(graph));
  
  Timer timer, totalTimer; timer.start(); totalTimer.start();

  auto outputPath = outputDir(pos.outputDir, false);
  FLTM_Data fltm_data;  
  std::cout << "Loading geno data from " << pos.inputDataFile << std::endl; // todo: logging
  fltm_data.matrix = loadDataTable(pos.inputDataFile);
  std::cout << "Loading pheno data from " << pos.inputPheno<< std::endl; // todo: logging
  auto pheno = loadPhenotype( pos.inputPheno );
  std::cout << "Loading label data from " << pos.inputLabelFile << std::endl; // todo: logging
  loadLabelPosition2( fltm_data.labels, fltm_data.indexes, fltm_data.positions, pos.inputLabelFile );
  
  std::cout << "Loading mapping snp data from " << pos.mappingFile << std::endl; // todo: logging
  std::map<std::string, std::string> snp2rs;
  snpMapping(snp2rs, pos.mappingFile);

  printf("assuring the graph position\n\n");
  assureGraphPositions(graph);
  std::cout << "done pre-processing takes: " << timer.display() << std::endl << std::endl;
  timer.restart();

  std::cout << "retrieving the candidates...\n\n";
  auto candidatesByLevel = retrieveCandidates(graph);

  std::cout << "retrieving the cardinalities...\n\n";
  auto cardinalities = retrieveCardinalities(graph);

  // auto outFileName = ( outputPath / "gwas_p_values_result.txt").string();
  auto chisq = std::make_shared<ChiSq>();

  std::cout << "start performing test...\n";

  std::vector<double> thresholds {
    0.5,1,1, 1, 1, 1, 1, 1
  };

  if ( pos.task == 0 ) { // gwas
    performHasGoodParentPermutationTest
        ( fltm_data, graph,
          *pheno, chisq,
          candidatesByLevel,
          cardinalities,
          pos.permutations,
          pos.chromosome,
          thresholds,
          outputPath );
  }
  else if (pos.task == 1) {
    performTest
        ( fltm_data, graph, snp2rs,
          *pheno, chisq,
          candidatesByLevel,
          cardinalities,
          pos.permutations,
          pos.chromosome,
          thresholds,
          outputPath );  
  } else {

  }
    // what we do here is:
    // we inforce that children should have parent in uppler
    // performTest
    //     ( fltm_data, graph,
    //              *pheno, chisq,
    //              candidatesByLevel,
    //              cardinalities,
    //              pos.permutations,
    //              pos.chromosome,
    //              pos.threshold,
    //              outputPath );
  
  // } else {
  //   // auto outFileName = ( outputPath / "regions.txt").string();
  //   // std::cout << "performing saving regions task to: " << outFileName << std::endl;
  //   // saveRegions( outFileName, pos.chromosome, graph );   
  // }
  std::cout << "done! program now exits\n";
}

//////////////////////////////////////////////////////////////
CandidatesByLevel retrieveCandidates( const Graph& graph ) {
  int max_level = -1;
         
  for ( auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first ) {
    auto vertex = *vp.first;
    auto &node = graph[vertex];
    max_level = std::max( node.level, max_level );
  }

  CandidatesByLevel candidates(max_level+1, Candidates());  
  for ( auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first ) {
    auto vertex = *vp.first;
    auto &node = graph[vertex];
    int level = node.level;
    candidates[level].push_back(vertex);
  }

  return candidates;
}

//////////////////////////////////////////////////////////////////////////

Cardinalities retrieveCardinalities( const Graph& graph ) {
  Cardinalities cards;
  for ( auto vp = boost::vertices(graph); vp.first != vp.second; ++vp.first ) {
    auto vertex = *vp.first;
    auto &node = graph[vertex];
    cards.push_back(node.variable.cardinality());
  }
  return cards;
}

/////////////////////////////////////////////////////////////////////////
Candidates getChildrenFromCandidates( const Graph& graph,
                                      const Candidates& curr_candidates ) {
  Candidates rs;
  for ( auto cand: curr_candidates ) {
    for ( auto ei = boost::out_edges(cand,graph); ei.first != ei.second; ++ei.first ) {
      auto child = boost::target(*ei.first, graph);
      rs.push_back(child);
    }  
  } 
  return rs;
}


/////////////////////////////////////////////////////////////////////////

void performTest( const FLTM_Data& fltmData, const Graph& graph, std::map<std::string, std::string>& snp2rs,
                  const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                  const CandidatesByLevel& candidatesByLevel,
                  const Cardinalities &cardinalities,
                  const int permutations,
                  const int chr,
                  const std::vector<double> thresholds,
                  boost::filesystem::path outDir) {

  auto gwas_outFile = ( outDir / "gwas.txt").string();
  auto parent = getGraphParent(graph); 

  char gwas_filtered_buffer[80], region_filtered_buffer[80];
  sprintf( gwas_filtered_buffer, "gwas_filtered.txt");
  sprintf( region_filtered_buffer, "region_filtered.txt");

  auto gwas_filtered_outFile = ( outDir / gwas_filtered_buffer).string();
  auto region_outFile = ( outDir / region_filtered_buffer).string();
  
  std::ofstream gwasFile(gwas_outFile), gwasFilteredFile(gwas_filtered_outFile), regionFile(region_outFile);
  printf("beginning of testing of %d permutations by test: %s\n", permutations, statTest->name.c_str());
  int levels = candidatesByLevel.size();

  std::vector<double> scores(fltmData.matrix->size(), 0.0);

  gwasFile << "chr" << chr << SEPARATOR
           << "id" << SEPARATOR
           << "snp" << SEPARATOR
           <<  "snp-id" << SEPARATOR
           << "level" << SEPARATOR
           << "parent" << SEPARATOR
           << "position"<< SEPARATOR
           << "pval" << SEPARATOR
           << "corrected-pval" << std::endl;
  
  for ( int l = 0; l < candidatesByLevel.size(); ++l) {
    auto &candidates = candidatesByLevel[l];
    auto cards = retrieveCardinalities(graph,candidates);
    
    std::vector<double> dist;
    std::vector<double> pvals;
    if ( candidates.size() ) {
      printf("\nWe now tests %d vars - @level: %d over %d\n\n", candidates.size(), l, levels - 1);
      permutationProcedure( dist, pvals, statTest,
                            *fltmData.matrix, pheno,
                            candidates, cards,
                            2, permutations);



      for ( int i = 0; i < candidates.size(); ++i ) {        
        auto cand = candidates[i];
        std::string rs = "";
        if ( snp2rs.find(graph[cand].label) != snp2rs.end() ) {
          rs =  snp2rs[graph[cand].label];
        } else {
          rs = "NA";
        }        
        scores[cand] =  pvals[2*i+1];
        gwasFile << "chr" << chr << SEPARATOR
                 << fltmData.indexes[cand] << SEPARATOR
                 << graph[cand].label << SEPARATOR
                 << rs << SEPARATOR
                 << l << SEPARATOR
                 << parent[cand] << SEPARATOR
                 << graph[cand].position << SEPARATOR
                 << pvals[2*i] << SEPARATOR
                 << pvals[2*i+1] << std::endl;

        if ( scores[cand] < thresholds[l]) {
          Node node = graph[cand];
          size_t sz_start = std::numeric_limits<int>::max(), sz_end = 0;
          size_t sum = 0, pos = 0, count = 0;
          for ( auto ei = boost::out_edges(cand, graph); ei.first != ei.second; ++ei.first ) {
            auto child = graph[boost::target(*ei.first, graph)];
            sz_start = std::min(sz_start, (size_t)child.position);
            sz_end = std::max(sz_end, (size_t)child.position);
            count++;
          }         
         
          gwasFilteredFile << "chr" << chr << SEPARATOR
                           << cand << SEPARATOR
                           << graph[cand].label << SEPARATOR
                           << l << SEPARATOR
                           << parent[cand] << SEPARATOR
                           << graph[cand].position << SEPARATOR
                           << pvals[2*i] << SEPARATOR
                           << pvals[2*i+1] << std::endl;

            
          if ( sz_start != sz_end && sz_start != node.position && sz_end != node.position  && count > 1 ) {
            printf("level: %d, var: %d, score: %f, thres: %f, start: %d, end: %d\n", l, cand, scores[cand], thresholds[l], sz_start, sz_end);
          if (sz_start < 0) {
            printf("level: %d, var: %d, score: %f, thres: %f, start: %d, end: %d\n", l, cand, scores[cand], thresholds[l], sz_start, sz_end);
            exit(-1);

          }
            if ( l > 0) {
              regionFile << "chr" << chr << " "
                         << sz_start << " "
                         << sz_end << std::endl;
            }
          }
        }
        
      } // candidate
    }
  }

  gwasFile.close(); gwasFilteredFile.close(); regionFile.close();
  std::cout << "writing to: " << gwas_outFile << std::endl;
  std::cout << "writing to: " << gwas_filtered_outFile << std::endl;
  std::cout << "writing to: " << region_outFile << std::endl;
}

/////////////////////

void performFilteredPermutationTest( const FLTM_Data& fltmData, const Graph& graph,
                                     const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                                     const CandidatesByLevel& candidatesByLevel,
                                     const Cardinalities &cardinalities,
                                     const int permutations,
                                     const int chr,
                                     const std::vector<double> thresholds,
                                     boost::filesystem::path outDir) {

  // auto gwas_outFile = ( outDir / "gwas.txt").string();

  char gwas_filtered_buffer[80], region_filtered_buffer[80];
  sprintf( gwas_filtered_buffer, "gwas_filtered.txt");
  sprintf( region_filtered_buffer, "region_filtered.txt");

  auto gwas_filtered_outFile = ( outDir / gwas_filtered_buffer).string();
  auto region_outFile = ( outDir / region_filtered_buffer).string();
  
  std::ofstream gwasFilteredFile(gwas_filtered_outFile), regionFile(region_outFile);
  printf("beginning of testing of %d permutations by test: %s\n", permutations, statTest->name.c_str());
  int levels = candidatesByLevel.size();

  Candidates candidates = candidatesByLevel[candidatesByLevel.size()-1];
  std::vector<double> scores(fltmData.matrix->size(), 0.0);
  while ( !candidates.empty() ) {
    auto cards = retrieveCardinalities(graph,candidates);
    std::vector<double> dist;
    std::vector<double> pvals;
    printf("\nWe now tests %d vars - @level over %d\n\n", candidates.size(), levels - 1);
    permutationProcedure( dist, pvals, statTest, *fltmData.matrix, pheno,
                          candidates, cards, 2, permutations);

    Candidates filtered_cands;
    for ( int i = 0; i < candidates.size(); ++i ) {
      auto cand = candidates[i]; auto &node = graph[cand]; auto l = node.level;
      scores[cand] =  pvals[2*i];
      double corrected = pvals[2*i+1];
      if (corrected<1)
        
      if ( corrected < thresholds[l] ) {
        filtered_cands.push_back(cand);
        Node node = graph[cand];
        int sz_start = std::numeric_limits<int>::max(), sz_end = -std::numeric_limits<int>::max();
        int sum = 0, pos = 0, count = 0;
        for ( auto ei = boost::out_edges(cand, graph); ei.first != ei.second; ++ei.first ) {
          auto child = graph[boost::target(*ei.first, graph)];
          sz_start = std::min(sz_start, child.position); sz_end = std::max(sz_end, child.position); count++;
        }          
        gwasFilteredFile << "chr" << chr << SEPARATOR << cand << SEPARATOR
                         << graph[cand].label << SEPARATOR << l << SEPARATOR
                         << graph[cand].position << SEPARATOR << pvals[2*i] << SEPARATOR
                         << pvals[2*i+1] << std::endl;
            
        if ( sz_start != sz_end && sz_start != node.position && sz_end != node.position  && count > 1 ) {
          printf("level: %d, var: %d, score: %f, start: %d, end: %d\n",
                 l, cand, scores[cand], sz_start, sz_end);
          if ( l > 0) {
            regionFile << "chr" << chr << " " << sz_start << " "   << sz_end << std::endl;
          }
        }
      }        
    }
    candidates = getChildrenFromCandidates(graph,filtered_cands);
  }

  gwasFilteredFile.close(); regionFile.close();
  // std::cout << "writing to: " << gwas_outFile << std::endl;
  std::cout << "writing to: " << gwas_filtered_outFile << std::endl;
  std::cout << "writing to: " << region_outFile << std::endl;
}




//////////////////p


Cardinalities retrieveCardinalities( const Graph& g, const Candidates& candidates ) {
  Cardinalities cards;
  for ( auto cand: candidates ) {
    cards.push_back(g[cand].variable.cardinality());
  }
  return cards;
}

///////////////////
bool verifyVariable( const Vertex v, const Graph& graph, const std::vector<int>& parent, const std::vector<double>& scores, double threshold ) {
  if ( scores[v] >= threshold ) return false;  
  auto latent = parent[v];
  if ( latent == -1 ) return true;
  return scores[latent] < threshold;
  // return scores[v] < threshold;
}

void performHasGoodParentPermutationTest( const FLTM_Data& fltmData, const Graph& graph,
                                          const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                                          const CandidatesByLevel& candidatesByLevel,
                                          const Cardinalities &cardinalities,
                                          const int permutations,
                                          const int chr,
                                          const std::vector<double> thresholds,
                                          boost::filesystem::path outDir) {
    // auto gwas_outFile = ( outDir / "gwas.txt").string();
  auto gwas_outFile = ( outDir / "gwas.txt").string(); std::ofstream gwasFile(gwas_outFile);
  
  char gwas_filtered_buffer[80], region_filtered_buffer[80];
  sprintf( gwas_filtered_buffer, "gwas_ok_parent_filtered.txt");
  sprintf( region_filtered_buffer, "region_ok_parent_filtered.txt");
  auto parent = getGraphParent(graph); 
  auto gwas_filtered_outFile = ( outDir / gwas_filtered_buffer).string();
  auto region_outFile = ( outDir / region_filtered_buffer).string();
  
  std::ofstream gwasFilteredFile(gwas_filtered_outFile), regionFile(region_outFile);
  printf("beginning of testing of %d permutations by test: %s\n", permutations, statTest->name.c_str());
  int levels = candidatesByLevel.size();

  std::vector<double> scores(fltmData.matrix->size(), 0.0);
  for ( int l = levels - 1; l >= 0; --l ) {
    auto &candidates = candidatesByLevel[l];
    auto cards = retrieveCardinalities(graph,candidates);
    std::vector<double> dist;
    std::vector<double> pvals;
    printf("\nWe now tests %d vars - @level %d over %d\n\n", candidates.size(), l, levels - 1);
    permutationProcedure( dist, pvals, statTest, *fltmData.matrix, pheno,
                          candidates, cards, 2, permutations);
    printf("done testing, now processing\n");
    for ( int i = 0; i < candidates.size(); ++i ) {
      auto cand = candidates[i]; auto &node = graph[cand]; auto l = node.level;
      scores[cand] =  pvals[2*i+1];
      double corrected = pvals[2*i+1];
      gwasFile << "chr" << chr << SEPARATOR << cand << SEPARATOR
                         << graph[cand].label << SEPARATOR << l << SEPARATOR
                         << parent[cand] << SEPARATOR        
                         << graph[cand].position << SEPARATOR << pvals[2*i] << SEPARATOR
                         << pvals[2*i+1] << std::endl;
          
      if ( verifyVariable(cand, graph, parent, scores, thresholds[l] ) ) {
        Node node = graph[cand];
        int sz_start = std::numeric_limits<int>::max(), sz_end = -std::numeric_limits<int>::max();
        int sum = 0, pos = 0, count = 0;
        for ( auto ei = boost::out_edges(cand, graph); ei.first != ei.second; ++ei.first ) {
          auto child = graph[boost::target(*ei.first, graph)];
          sz_start = std::min(sz_start, child.position); sz_end = std::max(sz_end, child.position); count++;
        }
        std::cout << "chr" << chr << SEPARATOR << cand << SEPARATOR
                         << graph[cand].label << SEPARATOR << l << SEPARATOR
                         << parent[cand] << SEPARATOR        
                         << graph[cand].position << SEPARATOR << pvals[2*i] << SEPARATOR
                         << pvals[2*i+1] << std::endl;
        gwasFilteredFile << "chr" << chr << SEPARATOR << cand << SEPARATOR
                         << graph[cand].label << SEPARATOR << l << SEPARATOR
                         << parent[cand] << SEPARATOR        
                         << graph[cand].position << SEPARATOR << pvals[2*i] << SEPARATOR
                         << pvals[2*i+1] << std::endl;
            
        if ( sz_start != sz_end && sz_start != node.position && sz_end != node.position  && count > 1 ) {
          printf("level: %d, var: %d, score: %f - %f, start: %d, end: %d\n",
                 l, cand, pvals[2*i], corrected, sz_start, sz_end);
          if ( l > 0) {
            regionFile << "chr" << chr << " " << sz_start << " "   << sz_end << std::endl;
          }
        }
      } // verify        
    } // candidates iter
  } // level

  gwasFilteredFile.close(); regionFile.close(); gwasFile.close();
  std::cout << "writing to: " << gwas_outFile << std::endl;  
  std::cout << "writing to: " << gwas_filtered_outFile << std::endl;
  std::cout << "writing to: " << region_outFile << std::endl;
}

void snpMapping(std::map<std::string, std::string>& snp2rs, std::string& infile) {

  std::ifstream mapFile(infile.c_str());
  if (!mapFile) {
    printf("file mapping %s not existing\n", infile.c_str());
    exit(-1);
  }
  
  utility::CSVIterator<std::string> mapLine(mapFile);// ++labPosLine;
  for( ; mapLine != utility::CSVIterator<std::string>(); ++mapLine ) {
    auto snp = (*mapLine)[0];
    auto rs = (*mapLine)[1];
    snp2rs[snp] = rs;
  }

  mapFile.close();
}

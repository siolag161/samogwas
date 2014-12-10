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

void performTest( const FLTM_Data& fltmData, const Graph& g,
                  const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                  const CandidatesByLevel& candidates,
                  const Cardinalities &cardinalities,
                  std::string& outFile );

/*  permutationProcedure( std::vector<double> &distri,
    std::vector<double> &pvals,
    StatTest* statTest,
    const Matrix &xData,
    const Vector &yData,
    const std::vector<size_t> &xCandidates,
    const std::vector<int> &xCardinalities,
    const int yCardinality,
    const int permu = 2)*/

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

  auto outputPath = outputDir(pos.outputDir);
  if ( pos.task == 0 ) { // gwas
    FLTM_Data fltm_data;  
    std::cout << "Loading geno data from " << pos.inputDataFile << std::endl; // todo: logging
    fltm_data.matrix = loadDataTable(pos.inputDataFile);
    std::cout << "Loading pheno data from " << pos.inputPheno<< std::endl; // todo: logging
    auto pheno = loadPhenotype( pos.inputPheno );
    std::cout << "Loading label data from " << pos.inputLabelFile << std::endl; // todo: logging
    loadLabelPosition( fltm_data.labels, fltm_data.indexes, fltm_data.positions, pos.inputLabelFile );

    printf("assures the graph position\n");
    assureGraphPositions(graph);
    timer.restart();

    std::cout << "done pre-processing takes: " << timer.display() << std::endl;
    
    auto candidatesByLevel = retrieveCandidates(graph);
    auto cardinalities = retrieveCardinalities(graph);

    auto outFileName = ( outputPath / "gwas_p_values_result.txt").string();
    auto chisq = std::make_shared<ChiSq>();
    performTest( fltm_data, graph,
                 *pheno, chisq,
                 candidatesByLevel,
                 cardinalities,
                 outFileName );
  
  } else {
    // auto outFileName = ( outputPath / "regions.txt").string();
    // std::cout << "performing saving regions task to: " << outFileName << std::endl;
    // saveRegions( outFileName, pos.chromosome, graph );   
  }
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
  
  CandidatesByLevel candidates(max_level, Candidates());  
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

void performTest( const FLTM_Data& fltmData, const Graph& g,
                  const PhenoVec& pheno, std::shared_ptr<StatTest> statTest,
                  const CandidatesByLevel& candidates,
                  const Cardinalities &cardinalities,
                  std::string& outFile ) {
  
}

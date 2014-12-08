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
#include "data_load.hpp"
#include "app_options.hpp"

using namespace utility;
using namespace samogwas;
 
int main( int argc, char** argv ) {
  Options pos = getProgramOptions( argc, argv );

  FLTM_Data fltm_data;
  
  fltm_data.matrix = loadDataTable(pos.inputDataFile);
  loadLabelPosition( fltm_data.labels, fltm_data.indexes, fltm_data.positions, pos.inputLabelFile );

  Graph graph;
  BayesGraphLoad graphLoad;
  graphLoad( graph,
             pos.labelInFile,
             pos.graphBayesVertexInFile,
             pos.graphBayesDistInFile );
  
  std::cout << "Loading geno data from " << progOpt.genoInFile << std::endl; // todo: logging
  loadDataTable ( genoMat, progOpt.genoInFile );
  std::cout << "Loading pheno data from " << progOpt.phenoInFile << std::endl; // todo: logging
  loadPhenotype( pheno, progOpt.phenoInFile );
  std::cout << "Loading label data from " << progOpt.labelInFile << std::endl; // todo: logging

  
  
  auto graph = SingleGraphLoad()(pos.graphFile);
  auto bayes = 
  BayesGraphSave()( *result.graph, outBayesVertex, outBayesDist ); 
}


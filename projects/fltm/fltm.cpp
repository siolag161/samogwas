/****************************************************************************************
 * File: fltm.cpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/07/2014
 ***************************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> // to obtain the program's name

#include "utils/csv_parser.hpp"
#include "utils/matrix_utils.hpp"
#include "utils/timer_utils.hpp"
#include "utils/option_printer.hpp" //for printing

#include "distance/comparable.hpp"
#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "distance/dissimilarity.hpp"
#include "fltm/fltm.hpp"
#include "em/naive_bayes_em.hpp"

#include "application_options.hpp"
#include "data_load.hpp"

namespace utl = utility;
using namespace samogwas;

typedef std::vector< std::vector<int> > Matrix; // We consider here only vector of int is relevant
typedef samogwas::MutInfoDissimilarity<Matrix> MutInfoDiss; // 
typedef samogwas::MutInfoSimilarity<Matrix> MutInfoSimi;

void checkInputFiles( std::string& path, std::string filename );
void saveImputedData( std::string dataPath, std::string labposPath,
                      const FLTM_Data& input,                      
                      Matrix& mat,
                      const FLTM_Result& resultFLTM );

void printOptions( Options& opt ); // 
AlgoClusteringInterface* getAlgoClust( FLTM_Data& input, Options& opt );

/////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
  utl::Timer timer; timer.start();
  
  Options pos = getProgramOptions( argc, argv );
  std::cout << "loading data...\n";
  checkInputFiles ( pos.inputDataFile, "data"); checkInputFiles( pos.inputLabelFile, "label file");
  FLTM_Data fltm_data;
  loadDataTable ( fltm_data.matrix, pos.inputDataFile ); 
  loadLabelPosition( fltm_data.labels, fltm_data.indexes, fltm_data.positions, pos.inputLabelFile );
  Matrix tmpMatrix = fltm_data.matrix;
  printf("Data loaded. Matrix of nbrVariables: (%d,%d). Took %s\n", utility::nrows( fltm_data.matrix ), utility::ncols( fltm_data.matrix ), timer.display().c_str());
  timer.restart();

  std::cout << "Performing FLTM...\n";
  AlgoClusteringInterface* algoClust = getAlgoClust( fltm_data, pos );
  FLTM_Result result;
  LinearCardinality emLC(pos.fltm_alpha, pos.fltm_beta, pos.fltm_maxCard);
  NaiveBayesEM *multiEM = new NaiveBayesEM( pos.fltm_nbrRestarts, pos.fltm_imputeMode );
  FLTM fltm(algoClust, emLC, multiEM);

  BayesGraphSave bayesGraphSave;
  ///////////////
  boost::filesystem::path outputPath = boost::filesystem::absolute(pos.outputFile);
  std::string outBayesVertex, outBayesDist, outImpDat, outImpLab, outGraph;
  boost::filesystem::create_directories(outputPath);
  char bayesVertex_fn[256], bayesDist_fn[256], imputedDat_fn[256], imputedLab_fn[256], graph_fn[256];
  sprintf( bayesVertex_fn, "fltm_%s_bayes_vertex.out", algoClust->name() );
  sprintf( bayesDist_fn, "fltm_%s_bayes_dist.out", algoClust->name() );
  sprintf( imputedDat_fn, "fltm_%s_imputedData.out", algoClust->name() );
  sprintf( imputedLab_fn, "fltm_%s_imputedLab.out", algoClust->name() );
  sprintf( graph_fn, "fltm_%s_graph.out", algoClust->name() );

  outBayesVertex = (outputPath / bayesVertex_fn).string(),
      outBayesDist = (outputPath / bayesDist_fn).string(),
      outImpDat = (outputPath / imputedDat_fn).string(),
      outImpLab = (outputPath / imputedLab_fn).string(),
      outGraph = (outputPath / graph_fn).string();
  
  // SingleGraphSave()( fltm_data.graph, outGraph );
  // bayesGraphSave( fltm_data.graph, outBayesVertex, outBayesDist );  
  saveImputedData( outImpDat, outImpLab, fltm_data, tmpMatrix, result );

  std::cout << "BYE...\n" << std::endl;

}

/** checks if input exists and exists on giving the error message
 *
 */
void checkInputFiles( std::string& path, std::string filename ) {
  if ( !boost::filesystem::exists( path ) )
  {
    std::cout << "Can't find " << filename << " at " << path << "! Program will now close." << std::endl;
    exit(-1);
  }
}


void printOptions( Options& opt ) {
  printf("inDatFile: %s,\
inLabFile: %s, \
verbose: %d, \
outFile: %s, \
clustAlgo: %d, \
dbscan(%d,%.2f), \
cast(%.2f), \
cardFunc(%.2f,%.2f,%d), \
nbrRestart: %d\n \
Opts: card: %d, steps: %d, maxDist: %d, simiThres: %.2f, emThres: %.2f, latentVarQualityThres: %.2f\n",
         opt.inputDataFile.c_str(),
         opt.inputLabelFile.c_str(),
         opt.verbose,
         opt.outputFile.c_str(),
         opt.clustAlgo,
         opt.dbscan_minPts,
         opt.dbscan_eps,
         opt.cast_cast,
         opt.fltm_alpha, opt.fltm_beta, opt.fltm_maxCard,
         opt.fltm_nbrRestarts,
         opt.fltm_opts.cardinality,
         opt.fltm_opts.nbrSteps,
         opt.fltm_maxDist,
         opt.fltm_simiThres,
         opt.fltm_opts.emThres,
         opt.fltm_opts.latentVarQualityThres);

}

AlgoClusteringInterface* getAlgoClust( FLTM_Data& input, Options& opt ) {
   AlgoClusteringInterface* algo;
   if ( opt.clustAlgo == 0 ) { // DBSCAN
     MutInfoDiss* diss = new MutInfoDiss( input.matrix, input.positions, opt.fltm_maxDist, opt.fltm_simiThres );  
     algo = new DBSCAN<MutInfoDiss>( diss, opt.dbscan_minPts, opt.dbscan_eps );
     printf("DBSCAN(%d,%.2f)\n", opt.dbscan_minPts, opt.dbscan_eps );
   } else {
     MutInfoSimi* simi = new MutInfoSimi( input.matrix, input.positions, opt.fltm_maxDist, opt.fltm_simiThres );  
     algo = new CAST<MutInfoSimi>(simi,opt.cast_cast);
     printf("CAST(%.2f)\n", opt.cast_cast );
   }

   return algo;
}

void saveImputedData( std::string dataPath, std::string labposPath,
                      const FLTM_Data& input,
                      Matrix& mat,
                      const FLTM_Result& result ) {
  //std::cout << "saving imputed data...\n" << std::endl;
  // 
  ////////////////////////////////
  std::ofstream matOut(dataPath);
  printf("first: saving data of nbrVariables: %d\n", (int)(utility::nrows(mat)));

  for(size_t row = 0; row < utility::nrows(mat); row++) {
    for(size_t col = 0; col < utility::ncols(mat) - 1; col++) {
      matOut << mat[row][col] << ",";
    }
    if ( utility::ncols(mat) > 0) {
      matOut << mat[row][utility::ncols(mat)-1];
    }
    matOut << std::endl;
  }
  printf("then: saving imputed data of nbrVariables: %d\n", (int)(utility::nrows(result.imputedData)));

  for(size_t row = 0; row < utility::nrows(result.imputedData); row++) {
    for(size_t col = 0; col < utility::ncols(result.imputedData) - 1; col++) {
      matOut << result.imputedData[row][col] << ",";
    }
    if ( utility::ncols(result.imputedData) > 0) {
      matOut << result.imputedData[row][utility::ncols(result.imputedData)-1];
    }
    matOut << std::endl;
  }
  matOut.close();
  ////////////////////////////////////
  std::ofstream labPosOut(labposPath);
  vertex_iterator vi, vi_end;
  int latId = input.indexes[input.indexes.size()-1];
  // for ( boost::tie(vi, vi_end) = boost::vertices(input.graph); vi != vi_end; ++vi ) {
  //   vertex_t vertex = *vi;
  //   if ( vertex < input.indexes.size() )
  //     labPosOut << input.indexes[vertex] << "," << input.graph[vertex].label << "," << input.graph[vertex].position
  //               << "," << input.graph[vertex].variable.cardinality() << std::endl;
  //   else {
  //     labPosOut << ++latId << "," << "\"imputed-" + input.graph[vertex].label << "\"," << input.graph[vertex].position
  //               << "," << input.graph[vertex].variable.cardinality() << std::endl;
  //   }
  // }

  labPosOut.close();
}
 

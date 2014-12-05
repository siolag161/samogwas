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
#include "distance/dissimilarity.hpp"
#include "distance/similarity.hpp"
#include "fltm/fltm.hpp"

#include "utils/option_printer.hpp"
#include "utils/custom_option_desc.hpp"

#include "utils/logs_utils.hpp"
#include "data_load.hpp"
#include "application_options.hpp"

using namespace utility;
using namespace samogwas;

typedef std::vector< std::vector<int> > Matrix; // We consider here only vector of int is relevant
typedef samogwas::MutInfoDissimilarity<Matrix> MutInfoDiss; // 
typedef samogwas::MutInfoSimilarity<Matrix> MutInfoSimi;

void checkInputFiles( std::string& path, std::string filename );
void saveImputedData( std::string dataPath, std::string labposPath,
                      const FLTM_Data& input,                      
                      Matrix& mat,
                      const FLTM_Result& resultFLTM );

char* current_date();
boost::filesystem::path outputDir( Options& progOpt );

// void printOptions( Options& opt ); // 
std::shared_ptr<AlgoClusteringInterface> getAlgoClust( FLTM_Data& input, Options& opt );
typedef measure<std::chrono::seconds> mea; // measure execution time in seconds
 
int main( int argc, char** argv ) {
  Options pos = getProgramOptions( argc, argv );
  FLTM_Data fltm_data;
  // mea::log_execution( "load_data", "second", loadDataTable, fltm_data.matrix, pos.inputDataFile, ',','"' );
  fltm_data.matrix = loadDataTable(pos.inputDataFile);
      
  loadLabelPosition( fltm_data.labels, fltm_data.indexes, fltm_data.positions, pos.inputLabelFile );
  auto tmpMatrix = *fltm_data.matrix;
  printf("Data loaded. Matrix of nbrVariables: (%d,%d)n", utility::nrows( *fltm_data.matrix ), utility::ncols( *fltm_data.matrix ));

  std::cout << "Performing FLTM...\n";
  auto algoClust = getAlgoClust( fltm_data, pos );
  FLTM_Result result;
  LinearCardinality emLC(pos.fltm_alpha, pos.fltm_beta, pos.fltm_maxCard);
  auto multiEM = std::make_shared<NaiveBayesEM>( pos.fltm_nbrRestarts, pos.fltm_imputeMode );

  FLTM fltm(algoClust, emLC, multiEM);

  mea::execution( fltm, result, fltm_data, pos.fltm_opts  );
  
  auto outputPath = outputDir(pos);
  std::string outBayesVertex, outBayesDist, outImpDat, outImpLab, outGraph;
  boost::filesystem::create_directories(outputPath);
  char bayesVertex_fn[256], bayesDist_fn[256], imputedDat_fn[256], imputedLab_fn[256], graph_fn[256];
  sprintf(bayesVertex_fn, "fltm_%s_bayes.vertex", algoClust->name() );
  sprintf(bayesDist_fn, "fltm_%s_bayes.dist", algoClust->name() );
  sprintf(imputedDat_fn, "fltm_%s_imputedData.algo", algoClust->name() );
  sprintf(imputedLab_fn, "fltm_%s_imputedLab.lab", algoClust->name() );
  sprintf(graph_fn, "fltm_%s.graph", algoClust->name() );

  outBayesVertex = (outputPath / bayesVertex_fn).string(),
      outBayesDist = (outputPath / bayesDist_fn).string(),
      outImpDat = (outputPath / imputedDat_fn).string(),
      outImpLab = (outputPath / imputedLab_fn).string(),
      outGraph = (outputPath / graph_fn).string();
  
  SingleGraphSave()( *result.graph, outGraph );
  BayesGraphSave()( *result.graph, outBayesVertex, outBayesDist );  
  saveImputedData( outImpDat, outImpLab, fltm_data, tmpMatrix, result );

  std::cout << "BYE...\n" << std::endl;

}


std::shared_ptr<AlgoClusteringInterface> getAlgoClust( FLTM_Data& input, Options& opt ) {
  std::shared_ptr<AlgoClusteringInterface> algo;
  if ( opt.clustAlgo == 0 ) { // DBSCAN
    auto diss = std::make_shared<MutInfoDiss>( input.matrix, input.positions, opt.fltm_maxDist, opt.fltm_simiThres );
    algo = std::make_shared<DBSCAN<MutInfoDiss>>( diss, opt.dbscan_minPts, opt.dbscan_eps );
    printf("DBSCAN(%d,%.2f)\n", opt.dbscan_minPts, opt.dbscan_eps );
  } else {
    auto simi = std::make_shared<MutInfoSimi>( input.matrix, input.positions, opt.fltm_maxDist, opt.fltm_simiThres );
    algo = std::make_shared<CAST<MutInfoSimi>>( simi, opt.cast_cast);
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

  auto& imputedData = *result.imputedData;
  printf("then: saving imputed data of nbrVariables: %d\n", (int)(utility::nrows(imputedData)));
  for(size_t row = 0; row < utility::nrows(imputedData); row++) {
    for(size_t col = 0; col < utility::ncols(imputedData) - 1; col++) {
      matOut << (imputedData)[row][col] << ",";
    }
    if ( utility::ncols(imputedData) > 0) {
      matOut << (imputedData)[row][utility::ncols(imputedData)-1];
    }
    matOut << std::endl;
  } 
  matOut.close();
  ////////////////////////////////////
  auto& graph = *result.graph;
  std::ofstream labPosOut(labposPath);
  vertex_iterator vi, vi_end;
  int latId = input.indexes[input.indexes.size()-1];
  for ( boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi ) {
    vertex_t vertex = *vi;
    if ( vertex < input.indexes.size() )
      labPosOut << input.indexes[vertex] << "," << graph[vertex].label << "," << graph[vertex].position
                << "," << graph[vertex].variable.cardinality() << std::endl;
    else {
      labPosOut << ++latId << "," << "\"imputed-" + graph[vertex].label << "\"," << graph[vertex].position
                << "," << graph[vertex].variable.cardinality() << std::endl;
    }
  }

  labPosOut.close();
}
 

char* current_date()
{
    time_t rawtime;
    struct tm * timeinfo;
    char* buffer = new char[80];
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (buffer,80,"%Y_%m_%d_%H_%M_%S",timeinfo);

    return buffer;
}

boost::filesystem::path outputDir( Options& progOpt ) {
  auto path = boost::filesystem::absolute(progOpt.outputDir);
  // boost::filesystem::create_directories(path);
  path /= current_date();
  boost::filesystem::create_directories(path);

  return path;
}

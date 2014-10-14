#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <map>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <math.h>
#include <boost/lexical_cast.hpp>

#include "em/naive_bayes_em.hpp"
#include "utils/csv_parser.hpp"
#include "utils/matrix_utils.hpp"
#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "distance/dissimilarity.hpp"
#include "fltm/fltm.hpp"

#include "data_generation.hpp"
using namespace samogwas;
using namespace utility;
using namespace data_gen;
class Data 
{ 
};
BOOST_FIXTURE_TEST_SUITE( Test_FLTM, Data ) 


void initOptions( FLTM_Options& opt ) {
  opt.cardinality = 3;
  opt.nbrSteps = 2;
  opt.emThres = 0.0000001;
  opt.latentVarQualityThres = 0.0;
}

void initData( FLTM_Data& input, int nclusts, int N, int ncols, int CARD) {
  int nrows = nclusts*N;
  for ( int i = 0; i < nrows; ++i ) {
    input.positions.push_back(i);
    input.labels.push_back(boost::lexical_cast<std::string>(i));
  }
  input.matrix = GenerateClusteredData( nclusts, N, CARD, ncols )();  
}

//////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( Test_Def ) {
  // define a clustering
  // plug in the fltm
  size_t nclusts = 5, ncols = 40;
  size_t N = 3, CARD = 3, MAX_POS = 50;
  int nrows = nclusts*N;

  FLTM_Data input;
  FLTM_Options opt;
  FLTM_Result result;
  initData( input, nclusts, N, ncols, CARD );
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  // auto data = GenerateClusteredData( nclusts, N, CARD, ncols )();
  typedef samogwas::MutInfoDissimilarity<data_gen::Matrix> MutInfoDiss;

  MutInfoDiss* diss = new MutInfoDiss( input.matrix, input.positions, MAX_POS, -1);
  
  initOptions(opt);

  LinearCardinality emLC(0.2, 1, 5);
  EMInterface * multiEM = new NaiveBayesEM(CARD, 3);

  AlgoClusteringInterface* dbscan = new DBSCAN<MutInfoDiss>( diss, 2, 0.2);
  FLTM fltm( dbscan, emLC, multiEM ); 
  fltm( result, input, opt );
  
}


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here






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

#include "test_clustering.hpp"

class Data 
{ 
};

BOOST_FIXTURE_TEST_SUITE( Test_Euclidian, Data ) 

Matrix generateEuclidianClusteredData( const size_t, const size_t, const int );

BOOST_AUTO_TEST_CASE( Test_DBSCAN_10 ) {
  
  size_t nclusts = 3, ncols = 20;
  size_t N = 5;
  int nrows = nclusts*N;
  auto data = generateEuclidianClusteredData( nclusts, N, ncols );

  for (int i = 0; i < nrows; ++i) {
    for (int j=0; j<ncols; ++j) {
      std::cout << data[i][j] << " ";
    }
    std::cout << std::endl;
  }
  // std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);

  auto* diss = new EucDiss(data);
  

  EucDBSCAN dbscan( diss, 2, 0.32 );
  Partition result = dbscan();
  for ( int i = 0; i < nrows; ++i ) {
     int expected_cluster = i / N;
     BOOST_CHECK_EQUAL(result.getLabel(i), expected_cluster );
  }

  printf("diff = %f\n", (*diss)(0,1));

}

Matrix generateEuclidianClusteredData( const size_t nclusts, const size_t N, const int ncols ) {
  size_t nbrVars = nclusts*N;
  Matrix result(nbrVars, std::vector<int>(ncols, 0.0));

  for ( size_t i = 0; i < nbrVars; ++i ) {
    // for ( size_t j = 0; j < ncols; ++ j ) {
    // result[i] = std::vector<int>(ncols, (i/N));
    // }
  }
  
  return result;
}


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here

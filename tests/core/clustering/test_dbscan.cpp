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

BOOST_FIXTURE_TEST_SUITE( Test_DBSCAN_gene, Data ) 

BOOST_AUTO_TEST_CASE( Test_DBSCAN_10 ) {
  size_t nclusts = 3, ncols = 20;
  size_t N = 5, CARD = 3, MAX_POS = 50;
  int nrows = nclusts*N;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = GenerateClusteredData( nclusts, N, CARD, ncols )();  
  MutInfoDiss* diss = new MutInfoDiss(data, positions, MAX_POS, -1);

  for ( int i = 0; i < nrows; ++i ) {
    for ( int j = i; j < nrows; ++j ) {
      // printf("diff(%d,%d)= %f\n", i, tao dajo naj, diss->compute( i,j ) );
    }
  }

  std::cout<< "-----------------" << std::endl;


  DBSCAN dbscan( diss, 2, 0.32 );
  Partition result = dbscan();
  for ( int i = 0; i < nrows; ++i ) {
     int expected_cluster = i / N;
       BOOST_CHECK_EQUAL(result.getLabel(i), expected_cluster );
  }

  printf("diff = %f\n", (*diss)(0,1));
}





BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here

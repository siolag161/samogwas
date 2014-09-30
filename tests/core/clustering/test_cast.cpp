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

BOOST_FIXTURE_TEST_SUITE( Test_CAST_gene, Data ) 

BOOST_AUTO_TEST_CASE( Test_Gaussian_K_CAST_10 ) {
  
  size_t nclusts = 5, ncols = 40;
  size_t N = 3, CARD = 3, MAX_POS = 50;
  int nrows = nclusts*N;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = GenerateClusteredData( nclusts, N, CARD, ncols )();  
  MutInfoSimi* diss = new MutInfoSimi(data, positions, MAX_POS, -1);  
  CAST cast( diss, 0.5 );  
  Partition result = cast();
  for ( int i = 0; i < nrows; ++i ) {
    int expected_cluster = i / 3;
    BOOST_CHECK_EQUAL(result.getLabel(i), expected_cluster );
  }
  
}





BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here

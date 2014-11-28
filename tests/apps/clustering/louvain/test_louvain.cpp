#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>

#include <omp.h>

#include "clustering/louvain/graph.hpp"
#include "clustering/louvain/louv.hpp"

#include "data_generation.hpp"
class Data 
{ 

};

using namespace samogwas;
using namespace samogwas::louvain;

BOOST_FIXTURE_TEST_SUITE( Test_Louvain, Data ) 

Graph generateGraph(int nbrNodes, int nbrComm);


BOOST_AUTO_TEST_CASE( Test_Dummy ) {

}

Graph generateGraph(int nbrNodes, int nbrComm) {
  // Graph g;

  //  size_t nclusts = 5, ncols = 40;
  // size_t N = 3, CARD = 3, MAX_POS = 50;
  // int nrows = nclusts*N;
  // std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  // auto data = GenerateClusteredData( nclusts, N, CARD, ncols )();  
  // MutInfoSimi* diss = new MutInfoSimi(data, positions, MAX_POS, -1);  
  // CAST cast( diss, 0.5 );  
  // Partition result = cast();
  // for ( int i = 0; i < nrows; ++i ) {
  //   int expected_cluster = i / 3;
  //   BOOST_CHECK_EQUAL(result.getLabel(i), expected_cluster );
  // }
  
  // return ;
}

BOOST_AUTO_TEST_SUITE_END()  /// 

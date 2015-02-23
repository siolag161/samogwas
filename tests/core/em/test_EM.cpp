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
#include "em/em.hpp"

#include "data_generation.hpp"
#include "clustering/cast.hpp"
#include "clustering/dbscan.hpp"
#include "distance/dissimilarity.hpp"
#include "fltm/graph.hpp"
#include "utils/matrix_utils.hpp"
class Data 
{ 
};

typedef samogwas::MutInfoDissimilarity<data_gen::Matrix> MutInfoDiss;
typedef samogwas::MutInfoSimilarity<data_gen::Matrix> MutInfoSimi;
typedef samogwas::DBSCAN<MutInfoDiss> DBSCAN;
typedef samogwas::CAST<MutInfoSimi> CAST;
typedef samogwas::Partition Partition;

using namespace data_gen;// using namespace samogwas;
using namespace utility;


BOOST_FIXTURE_TEST_SUITE( Test_EM, Data ) 

inline void prepareEM( MatrixPtr& emMat,
                       Variables& vars,
                       const MatrixPtr cltMat,
                       const samogwas::Graph& graph,
                       const std::vector<int>& cluster,
                       const std::vector<int> local2Global);



Variable createLatentVar( const int lab, const int cardinality ) {
  return Variable( boost::lexical_cast<std::string>(lab), plIntegerType(0, cardinality - 1) );
}

BOOST_AUTO_TEST_CASE( Test_EM_functional ) {
  
  size_t nclusts = 5, ncols = 40;
  size_t N = 3, CARD = 3, MAX_POS = 50;
  int nrows = nclusts*N;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = GenerateClusteredData( nclusts, N, CARD, ncols )();  
  // MutInfoSimi* diss = new MutInfoSimi(data, positions, MAX_POS, -1);
  auto diss = std::make_shared<MutInfoSimi>(data, positions, MAX_POS, -1);  

  CAST cast( diss, 0.5 );  
  Partition result = cast();
  
  for ( int i = 0; i < nrows; ++i ) {
    int expected_cluster = i / 3;
    BOOST_CHECK_EQUAL(result.getLabel(i), expected_cluster );
  }

  
  std::vector<int> local2Global;

  samogwas::Graph graph;
  for (size_t i = 0; i < N*nclusts; ++i) {
    samogwas::vertex_t vertex = createVertex( graph, 3, true,
                                    boost::lexical_cast<std::string>(i), i, 0 );

    local2Global.push_back(vertex) ;
  }

  for ( auto& clt: result.to_clustering() ) {
    auto emMat = std::make_shared<Matrix>();
    Variables vars;
    prepareEM( emMat, vars, data, graph, clt, local2Global );
    Variable latentVar = createLatentVar( boost::num_vertices(graph), 3 );
    samogwas::ResultEM resultEM;
    samogwas::NaiveBayesEM multiEM(3,1);
    multiEM( resultEM, latentVar, vars, emMat, 0.000001 );
  }
}


inline void prepareEM( MatrixPtr& emMat,
                Variables& vars,
                const MatrixPtr cltMat,
                const samogwas::Graph& graph,
                const std::vector<int>& cluster,
                const std::vector<int> local2Global ) {
  
  auto tEMMat = std::make_shared<Matrix>();
  tEMMat->reserve(cluster.size());
  vars.clear();  
  tEMMat->push_back( std::vector<int>( ncols(*cltMat), -1) );  
  for ( auto& it: cluster ) {    
    vars ^= graph[local2Global.at(it)].variable;
    tEMMat->push_back( cltMat->at(it) );      
  }
  emMat = Transpose(*tEMMat);
   // delete(tEMMat);
}


BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here

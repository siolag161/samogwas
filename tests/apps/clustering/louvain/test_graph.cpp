#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>

#include <omp.h>

#include "clustering/louvain/graph.hpp"
#include "clustering/louvain/community.hpp"
#include "clustering/louvain/louv.hpp"


#include "community_generation.hpp"
class Data 
{ 

};

using namespace samogwas;
using namespace samogwas::louvain;

BOOST_FIXTURE_TEST_SUITE( Test_Louvain_Graph, Data ) 

BOOST_AUTO_TEST_CASE( Test_Construction ) {
  int commCard = 5, commCount = 4;
  int nrows = commCard*commCount;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = data_gen::GenerateClusteredData( commCount, commCard, CARD, NCOLS )();  
  std::shared_ptr<SimilarityMatrix> simi(new MutInfoSimi(data, positions, MAX_POS, -1));

  Graph* g  =  new Graph(simi); 

  BOOST_CHECK_EQUAL( g->nbrNodes(), commCard*commCount );
  
}

std::vector<NodeIndex> adjNodes( Graph& g, const NodeIndex& i );

BOOST_AUTO_TEST_CASE( Test_Adjacent ) {
  int commCard = 5, commCount = 4;
  int nrows = commCard*commCount;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = data_gen::GenerateClusteredData( commCount, commCard, CARD, NCOLS )();  
  std::shared_ptr<SimilarityMatrix> simi(new MutInfoSimi(data, positions, MAX_POS, -1));
  Graph* g  =  new Graph(simi);

  int i = 3;
  auto adj_nodes =adjNodes(*g,i);
  BOOST_CHECK_EQUAL(adj_nodes.size(), i+MAX_POS-1);
  
  i = 2;
  adj_nodes =adjNodes(*g,i);
  BOOST_CHECK_EQUAL(adj_nodes.size(), i+MAX_POS-1);

  i = 7;
  adj_nodes =adjNodes(*g,i);
  BOOST_CHECK_EQUAL(adj_nodes.size(), 2*(MAX_POS-1));

}



std::vector<NodeIndex> adjNodes( Graph& g, const NodeIndex& i) {
  std::vector<NodeIndex> rs;
  Graph::AdjIte vi, ve;
  for ( boost::tie(vi,ve) = g.adjacentNodes(i); vi != ve; ++vi ) {
    // printf("%d-%d\n", i,*vi);
    rs.push_back(*vi);
  }
  return rs;
}
BOOST_AUTO_TEST_SUITE_END()  /// 

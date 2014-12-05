#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE
#endif  
#include <boost/test/unit_test.hpp>

#include <omp.h>

#include "clustering/louvain/graph.hpp"
#include "clustering/louvain/community.hpp"
// #include "clustering/louvain/louv.hpp"
#include "community_generation.hpp"

#include "clustering/louvain/louv.hpp"
#include "test_common.hpp"

class Data 
{ 

};

using namespace samogwas;
using namespace samogwas::louvain;

BOOST_FIXTURE_TEST_SUITE( Test_Louvain_Method, Data ) 


struct Simi: public SimilarityMatrix {
  Simi(std::vector< std::vector<double> > d): sim(d) {}
  virtual double compute( const size_t varA, const size_t varB ) { return sim[varA][varB]; }

  virtual size_t nbrVariables() const {
    return sim.size();
  }

  virtual void invalidate() {}

  std::vector< std::vector<double> > sim;
};

BOOST_AUTO_TEST_CASE( Test_Linked_Weights_Wikipedia ) {
  // http://en.wikipedia.org/wiki/Modularity_%28networks%29
  std::vector< std::vector<double> > sim {
    {0,1,1,0,0,0,0,0,0,1}, // 0
    {1,0,1,0,0,0,0,0,0,0}, // 1
    {1,1,0,0,0,0,0,0,0,0}, // 2 
    {0,0,0,0,1,1,0,0,0,1}, // 3
    {0,0,0,1,0,1,0,0,0,0}, // 4
    {0,0,0,1,1,0,0,0,0,0}, // 5
    {0,0,0,0,0,0,0,1,1,1},
    {0,0,0,0,0,0,1,0,1,0},
    {0,0,0,0,0,0,1,1,0,0},
    {1,0,0,1,0,0,1,0,0,0}    
  };
  
  std::shared_ptr<SimilarityMatrix> simi(new Simi(sim));
  auto louv = std::make_shared<MethodLouvain>(simi);

  louv->first_phase();
  double modul_1st = louv->network->modularity();
  louv->second_phase();
  double modul_2nd = louv->network->modularity();
  BOOST_CHECK_EQUAL(modul_1st , modul_2nd);
  auto clustering = louv->run();
  
  BOOST_CHECK_EQUAL( clustering.nbrClusters(), 3 );
  BOOST_CHECK_EQUAL( clustering.getLabel(0), clustering.getLabel(1));
  BOOST_CHECK_EQUAL( clustering.getLabel(1), clustering.getLabel(2));
  BOOST_CHECK_EQUAL( clustering.getLabel(2), clustering.getLabel(9));

  
  for ( NodeIndex n = 0; n < clustering.nbrItems(); ++n ) {
    if ( n < 3 ) BOOST_CHECK_EQUAL( clustering.getLabel(n), 0);
    else if ( n < 6 ) BOOST_CHECK_EQUAL( clustering.getLabel(n), 1);
    else if ( n < 9 ) BOOST_CHECK_EQUAL( clustering.getLabel(n), 2);
    else if ( n == 9 ) BOOST_CHECK_EQUAL( clustering.getLabel(n), 0);
  }
  
}


BOOST_AUTO_TEST_CASE(TEST_GENERATE) {
  int commCard = 5, commCount = 4;
  int nrows = commCard*commCount;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = data_gen::GenerateClusteredData( commCount, commCard, CARD, NCOLS )();  
  std::shared_ptr<SimilarityMatrix> simi( new MutInfoSimi(data, positions, MAX_POS, +1) );
  std::shared_ptr<Graph> g(new Graph(simi));


  std::ofstream wg("./input/wg_real.txt");
  auto sz = data->size();
  for ( int r = 0; r < sz; ++r) {
    for ( int c = r+1; c < sz; ++ c) {
      double w = simi->compute(c,r);
      if ( w > 0)
        wg << r << " " << c << " " << w << std::endl;
    }
  }
  
  Network network(g);
  
  BOOST_CHECK_EQUAL(network.nbrCommunities(), commCount*commCard); // initially
  BOOST_CHECK_EQUAL(network.nbrNodes(), commCount*commCard);
  BOOST_CHECK_EQUAL(network.modularity(), -0.5);

  for (int i = 0; i < commCard*commCount; ++i) {
    BOOST_CHECK_EQUAL( network.getCommunity(i), i ); // initially
    BOOST_CHECK_EQUAL( network.membersOf(i).size(), 1 ); // initially
    BOOST_CHECK_EQUAL( network.membersOf(i)[0], i ); // initially
  }

  auto louv = std::make_shared<MethodLouvain>(simi);
  auto clustering = louv->run();

  for ( NodeIndex n = 0; n < clustering.nbrItems(); ++n ) {
    BOOST_CHECK_EQUAL( clustering.getLabel(n), n/5 ); // initially
  }

  // 

}

BOOST_AUTO_TEST_SUITE_END()  /// 







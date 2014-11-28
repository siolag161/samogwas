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

#include "test_common.hpp"

class Data 
{ 

};

using namespace samogwas;
using namespace samogwas::louvain;

BOOST_FIXTURE_TEST_SUITE( Test_Louvain_Community, Data ) 

BOOST_AUTO_TEST_CASE( Test_Construction ) {
  int commCard = 5, commCount = 4;
  int nrows = commCard*commCount;
  std::vector<int> positions; for ( int i = 0; i < nrows; ++i ) positions.push_back(i);
  auto data = data_gen::GenerateClusteredData( commCount, commCard, CARD, NCOLS )();  
  std::shared_ptr<SimilarityMatrix> simi( new MutInfoSimi(data, positions, MAX_POS, -1) );
  std::shared_ptr<Graph> g(new Graph(simi));

  Network network(g);
  
  BOOST_CHECK_EQUAL(network.nbrCommunities(), commCount*commCard); // initially
  BOOST_CHECK_EQUAL(network.nbrNodes(), commCount*commCard);

  for (int i = 0; i < commCard*commCount; ++i) {
    BOOST_CHECK_EQUAL( network.communityOf(i), i ); // initially
    BOOST_CHECK_EQUAL( network.membersOf(i).size(), 1 ); // initially
    BOOST_CHECK_EQUAL( network.membersOf(i)[0], i ); // initially
  }

  // BOOST_CHECK_EQUAL( network.modularity(), 0.0 ); // initially

}

BOOST_AUTO_TEST_CASE( Test_Total_Modulariy_n ) {
  std::vector< std::vector<double> > sim {
    {0,0,0,0,0},
    {0,0,0,0,0},
    {0,0,0,0,0},
    {0,0,0,0,0}, 
    {0,0,0,0,0}       
  };
  
  // SimilarityMatrix* simi = new Simi(sim);
  std::shared_ptr<SimilarityMatrix> simi(new Simi(sim));

  std::shared_ptr<Graph> g(new Graph(simi));
  Network network(g);
  BOOST_CHECK_EQUAL( network.modularity(), -2.0);
}


BOOST_AUTO_TEST_CASE( Test_Linked_Weights ) {
  std::vector< std::vector<double> > sim {
    {0,1,2,0,0},
    {1,0,0,3,0},
    {2,0,0,3,4},
    {0,3,3,0,0}, 
    {0,0,4,0,0}       
  };
  
  // SimilarityMatrix* simi = new Simi(sim);
  std::shared_ptr<SimilarityMatrix> simi(new Simi(sim));

  std::shared_ptr<Graph> g(new Graph(simi));

  Network network(g);
  BOOST_CHECK_EQUAL( g->linkedWeights(0), 3.0);
  BOOST_CHECK_EQUAL( g->linkedWeights(1), 4.0);
  BOOST_CHECK_EQUAL( g->linkedWeights(2), 9.0);
  BOOST_CHECK_EQUAL( g->linkedWeights(3), 6.0);
  BOOST_CHECK_EQUAL( g->linkedWeights(4), 4.0);

  BOOST_CHECK_EQUAL( network.totalWeights(), 13);

  double total_w_r = 0.0;
  for ( int i = 0; i < network.nbrNodes(); ++i ) total_w_r += g->linkedWeights(i);
  BOOST_CHECK_EQUAL( 2*network.totalWeights(), total_w_r); // becasue 

  BOOST_CHECK( network.modularity() > -0.5  );

}

BOOST_AUTO_TEST_CASE( Test_Linked_Weights_Wikipedia ) {
  // http://en.wikipedia.org/wiki/Modularity_%28networks%29
  std::vector< std::vector<double> > sim {
    {0,1,1,0,0,0,0,0,0,1},
    {1,0,1,0,0,0,0,0,0,0},
    {1,1,0,0,0,0,0,0,0,0},
    {0,0,0,0,1,1,0,0,0,1},
    {0,0,0,1,0,1,0,0,0,0},
    {0,0,0,1,1,0,0,0,0,0},
    {0,0,0,0,0,0,0,1,1,1},
    {0,0,0,0,0,0,1,0,1,0},
    {0,0,0,0,0,0,1,1,0,0},
    {1,0,0,1,0,0,1,0,0,0}    
  };
  
  // SimilarityMatrix* simi = new Simi(sim);
  std::shared_ptr<SimilarityMatrix> simi(new Simi(sim));

  std::shared_ptr<Graph> g(new Graph(simi));
  Network network(g);

  // check graph proprties
  BOOST_CHECK_EQUAL( network.nbrNodes(), 10);
  BOOST_CHECK_EQUAL( network.nbrLinks(), 12);

  // BOOST_CHECK_EQUAL( g->linkedWeights(0) 120);
  for (int i=0;i<10;++i)
  {

    BOOST_CHECK_EQUAL( g->linkedWeights(i), network.tot_linked_weights[i]);
 
  }
  
  const int A = 0, B = 1, C = 2;

  network.moveNode(0, A);


  double b4 = network.modularity(), gain = network.modularityGain(1,A, network.sharedWeights(1,A)), loss = network.modularityLoss(1)/12;
  BOOST_CHECK_EQUAL( gain, 0.0625);
  network.moveNode(1, A);

  b4 = network.modularity(); gain = network.modularityGain(2,A)/12, loss = network.modularityLoss(2)/12;
  BOOST_CHECK_CLOSE( network.modularityGain(2,A, network.sharedWeights(2,A)), 0.131944445, 0.0001);
  network.moveNode(2, A);
  
  BOOST_CHECK_CLOSE( network.modularityGain(9,A), 0.01041667, 0.0001);
  network.moveNode(9, A);

  // BOOST_CHECK( network.modularityGain(9,A) > 0);
  BOOST_CHECK_CLOSE( network.modularityGain(3,B), 0.000, 0.000);
  network.moveNode(3, B);

  BOOST_CHECK_CLOSE( network.modularityGain(4,B), 0.0625, 0.000);
  network.moveNode(4, B);

  BOOST_CHECK_CLOSE( network.modularityGain(5,B), 0.1319444, 0.0001);
  network.moveNode(5, B);

  BOOST_CHECK_CLOSE( network.modularityGain(6,C), 0.0, 0.000);
  network.moveNode(6, C);

  BOOST_CHECK_CLOSE( network.modularityGain(7,C), 0.0625, 0.000);
  network.moveNode(7, C);

  BOOST_CHECK_CLOSE( network.modularityGain(8,C), 0.1319444, 0.0001);
  network.moveNode(8, C);

  BOOST_CHECK_CLOSE( network.modularityGain(8,A), -0.069444, 0.001);
  network.moveNode(8, A);
  network.moveNode(8, C);


  for (int i=0;i<10;++i)
  {
    if ( i < 3 || i == 9)
      BOOST_CHECK_EQUAL( network.communityOf(i), A);
    else if ( i < 6)
      BOOST_CHECK_EQUAL( network.communityOf(i), B);
    else if ( i < 9)
      BOOST_CHECK_EQUAL( network.communityOf(i), C);
  }


  
  BOOST_CHECK_EQUAL( network.nbrCommunities(), 3);
  BOOST_CHECK_EQUAL( network.totalWeights(), 12.0);
  BOOST_CHECK_CLOSE( network.modularity(), 0.4895833, 0.001);

  BOOST_CHECK_EQUAL( network.sharedWeights(0,B), 0.0 );
  BOOST_CHECK_EQUAL( network.sharedWeights(1,B), 0.0 );
  BOOST_CHECK_EQUAL( network.sharedWeights(9,B), 1.0 );
  BOOST_CHECK_EQUAL( network.sharedWeights(9,C), 1.0 );

  //BOOST_CHECK_EQUAL( network.modularityGain(9,C,-1.0) , 0.0 );
  BOOST_CHECK_EQUAL( network.tot_linked_weights[C] , 7 );
  BOOST_CHECK_EQUAL( network.linkedWeights(9), 3 );

  
  BOOST_CHECK_CLOSE( network.modularity(), 0.4895833, 0.001);


}

BOOST_AUTO_TEST_SUITE_END()  /// 

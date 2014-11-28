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
louv->run();
    // Network network(graph);

  // while (changed) {
  //   first_phase(network);
  //   second_phase(network);
  // }
  // auto graph = std::make_shared<Graph>(simi);
  // Network ntw(graph);
  // louv->first_phase(ntw);
  
}

BOOST_AUTO_TEST_SUITE_END()  /// 







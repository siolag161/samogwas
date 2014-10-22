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

#include "statistics/permutation_test.hpp"
#include "gwas/gwas_basic_strategy.hpp"
#include "gwas/score_node_criteria.hpp"

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <memory>
#include <boost/lockfree/queue.hpp>
using namespace stats;
using namespace samogwas;
using namespace std;
class Data
{ 
};

static Graph getGraph() {
  Graph g;
  for ( int i = 0; i <= 10; ++i ) boost::add_vertex(g);
  boost::add_edge(7,0,g); boost::add_edge(7,1,g); boost::add_edge(7,2,g);
  boost::add_edge(8,3,g); boost::add_edge(8,4,g);
  boost::add_edge(9,5,g); boost::add_edge(9,6,g);
  boost::add_edge(10,7,g); boost::add_edge(10,8,g); boost::add_edge(10,9,g); 
  return g;
}

BOOST_FIXTURE_TEST_SUITE( Test_Node_Criteria, Data )
//////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( Test_Threshold_Crtirion ) {
  // shared_ptr<ScoreCriterion> criterion_1(new LevelScoreNodeCriterion());

  Graph g = getGraph();
  
  std::vector<double> levelThres{ 0.05, 0.1, 0.1 };
  std::vector<double> scores{ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.01, 0.01};

  shared_ptr<NodeCriterion> criterion_1(new LevelScoreNodeCriterion(levelThres, scores));

  // BOOST_CHECK_EQUAL(criterion_1->isValid(graph,)
  for (int  i = 0; i < 11; ++i) {
    bool expected = ( i != 7 );
    BOOST_CHECK_EQUAL( expected, criterion_1->isValid(g,i) );
    expected = ( i == 7 );
    BOOST_CHECK_EQUAL( expected, criterion_1->isValid(g,i,false) );
  }
}


using namespace boost;

class bfs_visitor: public boost::default_bfs_visitor {
 public:  
  typedef std::map<Vertex, boost::default_color_type> ColorMap;
  typedef boost::associative_property_map<ColorMap> Color;
  bfs_visitor(std::shared_ptr<NodeCriterion> eval, Color& col): evaluator(eval), color(col) {}
  
  void initialize_vertex( Vertex v, const Graph &g) {
    std::cout << "Initialize: " << g[v].label << std::endl;
  }

  
  void examine_vertex( Vertex& v, const Graph &g ) {
    // std::cout << "exam: " << v << " " << color[v] <<  std::endl;
    printf("exam[%d]: color: %d\n", (int)v,(int)color[v]); 

  }
  
  void discover_vertex( Vertex& v, const Graph &g ) {
    printf("discover[%d]: color: %d\n", (int)v,(int)color[v]);
    typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      Vertex u = target(*ei, g);  
      color[u] = boost::color_traits<boost::default_color_type>::black();
    }
  }

 private:
  std::shared_ptr<NodeCriterion> evaluator;
  Color& color;
}; 



BOOST_AUTO_TEST_CASE( Test_Graph_Traversal ) {
  boost::queue<vertex_t> q; 
  typedef graph_traits<Graph>::vertex_descriptor vertex_t;
  typedef std::map<vertex_t, boost::default_color_type> ColorMap;
  typedef boost::associative_property_map<ColorMap> Color;
  Graph g = getGraph();
  ColorMap cmap; Color colorMap(cmap);

  
  std::vector<double> levelThres{ 0.05, 0.1, 0.1 };
  std::vector<double> scores{ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.01, 0.01};
  std::shared_ptr<NodeCriterion> criterion_1(new LevelScoreNodeCriterion(levelThres, scores));

  bfs_visitor vis(criterion_1, colorMap);//(&reached);
  boost::breadth_first_visit( g,
                              10,
                              q,
                              vis,
                              colorMap);
  // std::cout << "\npart deux: \n";

  // boost::queue<vertex_t> q1; 
  // ColorMap cmap1;
  // Color colorMap1(cmap1);
  // bfs_visitor vis1(criterion_1, colorMap);//(&reached);

  // boost::breadth_first_visit( g,
  //                             10,
  //                             q,
  //                             vis,
  //                             colorMap);


  // for ( int i = 0; i < 11; ++i ) {
  //   printf("color[%d] = %d\n", i, colorMap[i]);
  // }
  // std::cout << "part deux: \n";
}

BOOST_AUTO_TEST_SUITE_END()  /// Test InfoTheo ends here






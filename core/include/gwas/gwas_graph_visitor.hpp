/****************************************************************************************
 * File: gwas_graph_visitor.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 22/10/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_GWAS_GRAPH_VISITOR_HPP
#define SAMOGWAS_GWAS_GRAPH_VISITOR_HPP

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/lockfree/queue.hpp>

#include <fltm/graph.hpp>
#include <functional> // std::less
#include "node_criteria.hpp"
namespace samogwas
{

// typedef GWAS_Visitor {

// }
typedef boost::default_bfs_visitor GWAS_BFS_Visitor;

class GWAS_Basic_Visitor: public GWAS_BFS_Visitor {
  
 public:
  typedef NodeCriterion<double, std::less<double> > Criterion;
  typedef std::map<Vertex, boost::default_color_type> ColorMap;
  typedef boost::associative_property_map<ColorMap> Color;
  typedef std::map<Vertex, double> ScoreMap;
  
  GWAS_Basic_Visitor( std::shared_ptr<Criterion> eval,
                      std::shared_ptr<ScoreMap> scoreMap,
                      std::shared_ptr<Color> col)
      : evaluator(eval), scores(scoreMap), color(col) {
    visited_vertices = std::make_shared<std::vector<Vertex>>();
  }  

  GWAS_Basic_Visitor( std::shared_ptr<Criterion> eval,
                      std::shared_ptr<Color> col)
      : evaluator(eval), color(col) {
    scores = std::make_shared<ScoreMap>();
  }

  GWAS_Basic_Visitor( std::shared_ptr<Criterion> eval )
      : evaluator(eval) {
    scores = std::make_shared<ScoreMap>();
    color = std::make_shared<Color>();
  }
  
  void discover_vertex( Vertex& v, const Graph &g ) {
    std::cout << "visit: " << v << std::endl;
    bool rs = (scores && evaluator->isValid(g,v,*scores)) || !scores && evaluator->isValid(g,v);
    if (!rs) {
      typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
        Vertex u = target(*ei, g);
        if ( u != v ) {
          (*color)[u] = boost::color_traits<boost::default_color_type>::black();
        } // else {
        //   (*color)[u] = boost::color_traits<boost::default_color_type>::black();
        // }
      }
    } else {
      visited_vertices->push_back(v);
    }
  }

  std::vector<Vertex>& visitedVertices() const {
    return *visited_vertices;
  }
  
  // private:
  std::shared_ptr<Criterion> evaluator;
  std::shared_ptr<ScoreMap> scores;  
  std::shared_ptr<Color> color;

  std::shared_ptr<std::vector<Vertex>> visited_vertices;
}; 


} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_GWAS_GRAPH_VISITOR_HPP

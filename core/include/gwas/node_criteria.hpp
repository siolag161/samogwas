/****************************************************************************************
 * File: node_criteria.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 21/10/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_NODE_CRITERIA_HPP
#define SAMOGWAS_NODE_CRITERIA_HPP

#include <memory>
#include <vector>
// #include <functional.h>

#include <fltm/graph.hpp>

namespace samogwas
{

struct GraphNodeCriterion {
  virtual bool isValid( const Graph& g, const vertex_t& vertex ) const = 0;  
};

template<class T, class Func>
struct NodeCriterion: public GraphNodeCriterion {

  typedef typename std::map<Vertex, T> ScoreMap;
  
  virtual T nodeValue( const Graph& g, const vertex_t& vertex) const = 0;

  virtual T referenceValue( const Graph& g, const vertex_t& vertex) const = 0;
  
  virtual bool isValid( const Graph& g, const vertex_t& vertex ) const {
    T nodeVal = nodeValue(g,vertex);
    T referenceVal = referenceValue(g,vertex);
    
    return Func()(nodeVal, referenceVal);
  }

  virtual bool isValid( const Graph& g, const vertex_t& vertex, ScoreMap& scores) const {
    if ( scores.find(vertex) == scores.end() )
      scores[vertex] = nodeValue(g,vertex);

    T nodeVal = scores[vertex];
    T referenceVal = referenceValue(g,vertex);

    return Func()(nodeVal, referenceVal);
  }
  
};

template<class T, class Func>
struct NodeCriteria: public NodeCriterion<T,Func> {
  typedef NodeCriterion<T,Func> Criterion;
  typedef std::shared_ptr<Criterion> CriterionPtr;
  
  virtual bool isValid( const Graph& graph, const vertex_t& vertex ) const {
    for ( auto& criterion: criteria ) {
      if (!criterion->isValid(graph,vertex)) {
        return false; // gotta be all true
      }
    }
    return true; // ok if and only if all true
  }

 protected:
  std::vector<CriterionPtr> criteria;
};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_NODE_CRITERIA_HPP

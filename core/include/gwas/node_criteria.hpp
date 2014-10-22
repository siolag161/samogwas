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

struct NodeCriterion {
  // typedef std::less<double> LESS;
  // typedef std::greater<double> GREATER;
  // virtual double 
  virtual bool isValid( const Graph& graph, const vertex_t& vertex, bool less = true ) const = 0;   
};

struct NodeCriteria: public NodeCriterion {



  virtual bool isValid( const Graph& graph, const vertex_t& vertex, bool less = true ) const {
    for ( auto& criterion: criteria ) {
      if (!criterion->isValid(graph,vertex,less)) {
        return false; // gotta be all true
      }
    }
    return true; // ok if and only if all true
  }

 protected:
  std::vector< std::shared_ptr<NodeCriterion> > criteria;
};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_NODE_CRITERIA_HPP

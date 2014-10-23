/****************************************************************************************
 * File: gwas_basic_strategy.hpp
 * Description: This module provides @todo.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 01/08/2014

 ***************************************************************************************/

#ifndef SAMOGWAS_BASIC_GWAS_STRATEGY_HPP
#define SAMOGWAS_BASIC_GWAS_STRATEGY_HPP


#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/lockfree/queue.hpp>

#include <memory>
#include <vector>

#include <fltm/core_fltm.hpp>
#include <statistics/association_test.hpp>

#include "gwas_strategy.hpp"
#include "gwas_graph_visitor.hpp"
namespace samogwas
{

/** The visitor for traversing the FLTM graph
 *
 */

class GWAS_Basic_Strategy: public GWAS_Strategy {
  
 public:
  GWAS_Basic_Strategy( std::shared_ptr<GWAS_Basic_Visitor>  vis): visitor(vis) {}  

  GWAS_Basic_Strategy() {}  

  GWAS_Basic_Strategy& setVisitor( std::shared_ptr<GWAS_Basic_Visitor> vis ) {
    visitor = vis;
    return *this;
  }
  
  virtual void execute( Graph& graph,
                        Matrix& genotype,
                        Vector& phenotype );

 private:
  std::shared_ptr<GWAS_Basic_Visitor> visitor; 
};

///////////////////////////////////////////////////////////

class GWAS_Strategy_Builder {
  virtual std::shared_ptr<GWAS_Strategy> build( std::shared_ptr<GraphNodeCriterion> criteria ) = 0;
};


class GWAS_Basic_Strategy_Builder: public GWAS_Strategy_Builder {
  typedef NodeCriterion<double, std::less<double> > Criterion;

  virtual std::shared_ptr<GWAS_Strategy> build( std::shared_ptr<GraphNodeCriterion> nodeCriteria ) {
    std::shared_ptr<GWAS_Strategy> result;

    auto criteria = std::dynamic_pointer_cast<Criterion>(nodeCriteria);
    if (criteria) {
      auto vis = std::make_shared<GWAS_Basic_Visitor>(criteria);
      result =  std::make_shared<GWAS_Basic_Strategy>(vis);
    }

    return result;
  }

};

} // namespace samogwas ends here.

  /**************************************** IMPLEMENTATION BELOW THIS POINT ****************************************/


/****************************************************************************************/
#endif // SAMOGWAS_CORE_FLTM_HPP

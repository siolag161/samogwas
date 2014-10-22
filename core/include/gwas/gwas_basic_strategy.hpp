/****************************************************************************************
 * File: gwas_basic_strategy.hpp
 * Description: This module provides @todo.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 01/08/2014

 ***************************************************************************************/

#ifndef SAMOGWAS_BASIC_GWAS_STRATEGY_HPP
#define SAMOGWAS_BASIC_GWAS_STRATEGY_HPP

#include <vector>
#include <fltm/core_fltm.hpp>
#include <statistics/association_test.hpp>

#include "gwas_strategy.hpp"

namespace samogwas
{

/** The visitor for traversing the FLTM graph
 *
 */
class GWAS_Basic_Visitor;

class GWAS_Basic_Strategy: public GWAS_Strategy {
  
 public:
  GWAS_Basic_Strategy( std::vector<double>& thres ): thresholds(thres) {

  }

  // virtual std::map<int, std::vector<vertex_t> > 
      
  virtual void execute( FLTM_Result& result, Matrix& genotype, Vector& phenotype, stats::StatTest* statTest );

 private:
  std::vector<double>& thresholds;
};

} // namespace samogwas ends here.

/**************************************** IMPLEMENTATION BELOW THIS POINT ****************************************/


/****************************************************************************************/
#endif // SAMOGWAS_CORE_FLTM_HPP

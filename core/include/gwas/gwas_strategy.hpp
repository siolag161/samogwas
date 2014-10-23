/****************************************************************************************
 * File: gwas_strategy.hpp
 * Description: This module provides common interfaces for GWAS strategies.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 01/08/2014

 ***************************************************************************************/

#ifndef SAMOGWAS_GWAS_STRATEGY_HPP
#define SAMOGWAS_GWAS_STRATEGY_HPP

#include <vector>
#include <memory> // std::shared_ptr
#include "fltm/core_fltm.hpp"
#include "node_criteria.hpp"

namespace samogwas
{
class GWAS_Strategy {
  
 public:
  typedef std::vector< std::vector<vertex_t> > Level2Vertices;
  typedef std::map< int, std::vector<int> > Latent2Children;
  
  typedef std::vector<int> Vector;
  typedef std::vector<Vector> Matrix;

  virtual Level2Vertices levels2Vertices( const Graph& graph);
  virtual Latent2Children latent2Children(const Graph& graph);
  virtual void execute( Graph& graph,
                        Matrix& genotype,
                        Vector& phenotype ) = 0;
};


} // namespace samogwas ends here.

/**************************************** IMPLEMENTATION BELOW THIS POINT ****************************************/


/****************************************************************************************/
#endif // SAMOGWAS_GWAS_STRATEGY_HPP

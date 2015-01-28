/****************************************************************************************
 * File: compare_measure.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/01/2015

 ***************************************************************************************/
#ifndef SAMOGWAS_COMPARE_MEASURE_HPP
#define SAMOGWAS_COMPARE_MEASURE_HPP

#include "partition.hpp"
namespace samogwas
{


struct ComparisonMeasure {
  //  
  virtual double compute( const Clustering& c1, const Clustering& c2 ) const = 0;
  //
  virtual double compute( const Partition& p1, const Partition& p2 ) const  {
    Clustering c1 = p1.to_clustering();
    Clustering c2 = p2.to_clustering();
    return compute( c1, c2 );
  }  
  //  
  virtual double operator()( const Clustering& c1, const Clustering& c2 ) const { return compute(c1,c2); }
  //
  virtual double operator()( const Partition& p1, const Partition& p2 ) const { return compute(p1,p2); }

  virtual std::string name() const = 0;
};


} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_COMPARE_MEASURE_HPP

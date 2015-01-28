/****************************************************************************************
 * File: mirkin.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 09/01/2015

 ***************************************************************************************/
#ifndef SAMOGWAS_MIRKIN_HPP
#define SAMOGWAS_MIRKIN_HPP

#include "compare_measure.hpp"
#include "counter.hpp"
namespace samogwas
{

struct Mirkin: public ComparisonMeasure {
  
  virtual double compute( const Clustering& c1, const Clustering& c2 ) const {
    size_t c1_sum2 = 0;
    size_t n = 0;
    for (size_t i=0; i < c1.size(); i++) {
      c1_sum2 += c1[i].size() * c1[i].size();
      n += c1[i].size();
    }
    size_t c2_sum2 = 0;
    for (size_t i=0; i < c2.size(); i++) {
      c2_sum2 += c2[i].size() * c2[i].size();    
    }  
    size_t c1c2_sum2 = 0;
    for (size_t i=0; i < c1.size(); i++) {
      for (size_t j=0; j < c2.size(); j++) {
        size_t size;
        std::set_intersection( c1[i].begin(), c1[i].end(), 
                               c2[j].begin(), c2[j].end(),
                               counter(size) );
        c1c2_sum2 += size * size;
      }
    }

    return ( c1_sum2 + c2_sum2 - (2 * c1c2_sum2) ) / (double)(n*n);
  }

  //
  virtual double compute( const Partition& p1, const Partition& p2 ) const {
    Clustering c1 = p1.to_clustering();
    Clustering c2 = p2.to_clustering();
    return compute( c1, c2 );
  }
  
  virtual std::string name() const { return "Mirkin"; }

};

} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_MIRKIN_HPP

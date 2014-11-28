/****************************************************************************************
 * File: link.cpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 26/11/2014

 ***************************************************************************************/
#include <memory>
#include "distance/comparable.hpp"
#include "clustering/louvain/link.hpp"
namespace samogwas
{

namespace louvain {

WeightedLink::WeightedLink(WeightPtr w): weights(w) {}

WeightedLink::WeightedLink( const size_t sz) {
  //weights = std::make_shared<Weight>( sz, std::vector<double>(sz, 0.0) ;
  weights = std::make_shared<Weights>();
  weights->reserve(sz);
  for ( size_t i = 0; i < sz; ++i ) {
    weights->push_back( std::vector<double>(i+1, 0.0) );    
  }
}

size_t WeightedLink::nbrVariables() const {
  return weights->size();
}

double WeightedLink::compute( const size_t varA, const size_t varB ) {
  if (varA < varB) return compute(varB, varA);
  return (*weights)[varA][varB];
}

void WeightedLink::setWeight( const size_t varA, const size_t varB, const double val ) {
  if (varA < varB) return setWeight(varB, varA, val);
  // if ( varA == varB )
  //   printf("\nselfLoop of %d: %f\n\n", varA, val);
  (*weights)[varA][varB] = val;
}

}
} // namespace samogwasends here. samogwas

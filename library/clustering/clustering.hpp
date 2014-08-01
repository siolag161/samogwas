/****************************************************************************************
 * File: clustering.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_CLUSTERING_HPP
#define SAMOGWAS_CLUSTERING_HPP

#include "partition.hpp"
namespace samogwas
{

template<typename CompareMatrix>
class AlgoClust {
 public:
  virtual Partition operator()() = 0;
  // virtual * getComp() { return comp; }
  // virtual void setComp( CompMatrix* c ) { comp = c; }
  virtual char* name() const = 0;
  virtual void invalidCache() { comp->invalidCache(); }

  AlgoClust( CompareMatrix* c): comp(c) {}
  virtual ~AlgoClust() { delete comp; }
 protected:
  CompareMatrix* comp;
};


} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_CLUSTERING_HPP

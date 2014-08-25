/****************************************************************************************
 * File: clustering.hpp
 * Description: Common interfaces for Clustering algorithms
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_CLUSTERING_HPP
#define SAMOGWAS_CLUSTERING_HPP

#include "partition.hpp"
namespace samogwas
{

/** Any clustering algorithm is derived from this based class.
 *  This provide common methods like invalide caching and the name of the derived method
 */
class AlgoClustering {
 public:
  virtual Partition operator()() = 0;
  virtual char* name() const = 0;
  virtual void invalidCache() = 0;
};


/** The AlgoClust @Todo: Name changing is a specific case of AlgoClustering where the algorihm
 *  categorizes objects using a similarity/distance matrix.
 */
template<typename CompareMatrix>
class AlgoClust: public AlgoClustering {
 public:
  virtual Partition operator()() = 0;
  virtual char* name() const = 0;
  virtual void invalidCache() { comp->invalidCache(); }

  AlgoClust( CompareMatrix* c): comp(c) {}
  virtual ~AlgoClust() { delete comp; }
  
 protected:
  // Matrix of similarity/distance
  CompareMatrix* comp;
};


} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_CLUSTERING_HPP

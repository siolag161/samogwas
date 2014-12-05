/****************************************************************************************
 * File: clustering.hpp-------------------@todo: change the name to algo_clustering.hpp
 * Description: Common interfaces for clustering algorithms
 * @author: CS Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_ALGO_CLUSTERING_HPP 
#define SAMOGWAS_ALGO_CLUSTERING_HPP

#include <memory>

#include "partition.hpp"

namespace samogwas
{

/** Any clustering algorithm is derived from this base class.
 */
class AlgoClusteringInterface { 
 public:
  virtual Partition operator()() { return run(); }
  virtual Partition run() = 0;
  virtual char* name() const = 0;
  virtual void invalidate() = 0; //@todo: invalidateCache
};


/** The AlgoClust is a specific case of AlgoClusteringInterface where the algorihm
 *  categorizes objects using a similarity/dissimilarity matrix.
 *  The template parameter CompareMatrix has to be derived from 
 *  the samogwas::CompMatrix structure (comparable.hpp). @todo: change comparable.hpp and CompMatrix
 */
template<typename CompareMatrix> // CS CompareMatrix should be ComparisonMatrix
class AlgoClust: public AlgoClusteringInterface {
 public:
  virtual Partition run() = 0;
  virtual char* name() const = 0;
  virtual void invalidate() { compMatrix->invalidate(); }

  AlgoClust( std::shared_ptr<CompareMatrix> c): compMatrix(c) {}
  // virtual ~AlgoClust() { delete compMatrix; }
  
 protected:
  // Matrix of similarities or dissimilarities
  // CompareMatrix* compMatrix;
  std::shared_ptr<CompareMatrix> compMatrix;
};


} // namespace samogwas ends here.

/****************************************************************************************/
#endif // SAMOGWAS_ALGO_CLUSTERING_HPP

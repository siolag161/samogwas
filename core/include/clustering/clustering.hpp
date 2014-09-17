/****************************************************************************************
 * File: clustering.hpp
 * Description: Common interfaces for clustering algorithms
 * @author: CS siolag161 (thanh.phan@outlook.com)
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_CLUSTERING_HPP // CS Directive name assigment is not consistent throughout the project - 
                                // see partition.hpp where the directive FLTM_CLUSTERING_HPP is used 
                                // (it should be SAMOGWAS_PARTITION)
#define SAMOGWAS_CLUSTERING_HPP

#include "partition.hpp"
namespace samogwas
{

/** Any clustering algorithm is derived from this base class.
 *  This provides common methods like invalid caching and the name of the derived method CS INCOMPREHENSIBLE
 */
class AlgoClustering { // CS Is this level mandatory? Many redundancies with class AlgoClust.
 public:
  virtual Partition operator()() { return run(); }
  virtual Partition run() = 0;
  virtual char* name() const = 0;
  virtual void invalidCache() = 0;

 // static Partition to_partition( const std::vector<int>& labels );

};


/** The AlgoClust @Todo: Name changing is a specific case of AlgoClustering where the algorihm
 *  categorizes objects using a similarity/distance matrix.
 * CS
 */
template<typename CompareMatrix> // CS CompareMatrix should be ComparisonMatrix
class AlgoClust: public AlgoClustering {
 public:
  virtual Partition run() = 0;
  virtual char* name() const = 0;
  virtual void invalidCache() { comp->invalidCache(); }

  AlgoClust( CompareMatrix* c): comp(c) {}
  virtual ~AlgoClust() { delete comp; }
  
 protected:
  // Matrix of similarities or distances
  CompareMatrix* comp; // CS Bad identifier - Confusing with comparator - compMat should be used instead.
};


} // namespace samogwas ends here. samogwas // CS INCOMPREHENSIBLE

/****************************************************************************************/
#endif // SAMOGWAS_CLUSTERING_HPP

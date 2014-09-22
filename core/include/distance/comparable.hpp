/****************************************************************************************
 * File: comparable.hpp
 * Description: provides the common interface for CompareMatrix
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet 
 * @date: 30/07/2014
 ***************************************************************************************/
#ifndef SAMOGWAS_COMPARABLE_HPP
#define SAMOGWAS_COMPARABLE_HPP
 
#include "stddef.h" // size_t
namespace samogwas
{

/** Common Interface for dissimilarity matrix and similarity matrix.
 *  Usage: CompMatrix(a,b): returns either the similarity or dissimilarity
 *  between two objects indexed by a and b. @todo: change the interface name
 */
struct CompMatrix { 
  
  /** Computes and returns either the similarity or dissimilarity between two objects  
   * indexed by varA and varB in a given matrix. This is a functor facility, that in turn calls the compute method.
   */
  virtual double operator()( const size_t varA, const size_t varB ) {
    return compute(varA, varB);
  }

  /** Computes and returns either the similarity or dissimilarity between two objects  
   * indexed by varA and varB in a given matrix. 
   */
  virtual double compute( const size_t varA, const size_t varB ) = 0;

  /** Returns the number of total elements actually stored in the matrix (which is possibly sparse).
   */
  virtual size_t size() const = 0;

  /** Invalidates current caching if any. Caching may be helpful when the computation of the matrix
   * is performed on the go (example: mutual information as similarity and entropies stored in the cache).
   *  @TODO: invalid-> invalidate
   */
  virtual void invalidCache() = 0;

};

/////////////////////////////////////////////////////

/** Common Interface for dissimilarity matrix
 *  Usage: DissimilarityMatrix(a,b): returns the dissimilarity
 *  between two objects indexed by a and b.
 */
struct DissimilarityMatrix: public CompMatrix {};

/////////////////////////////////////////////////////

/** Common Interface for similarity matrix
 *  Usage: SimilarityMatrix(a,b): returns the similarity
 *  between two objects indexed by a and b.
 */
struct SimilarityMatrix: public CompMatrix {};

} // namespace samogwas ends here.

/****************************************************************************************/
#endif // SAMOGWAS_COMPARABLE_HPP

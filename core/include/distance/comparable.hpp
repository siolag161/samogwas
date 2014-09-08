/****************************************************************************************
 * File: comparable.hpp
 * Description: Provides the common interface for CompareMatrix
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_COMPARABLE_HPP
#define SAMOGWAS_COMPARABLE_HPP
 
#include "stddef.h" // size_t
namespace samogwas
{

/** Common Interface for Distance Matrix and Similarity Matrix.
 *  Usage: CompMatrix(a,b): return either the similarity or dissimilarity
 *  between 2 objects indexed by a and b
 */
struct CompMatrix {
  
  /** Computes and returns either the similarity or distance between 2 object
   * indexed by varA and varB in a given matrix. This is a proxy method, in turns call the compute method
   *
   */
  virtual double operator()( const size_t varA, const size_t varB ) {
    return compute(varA, varB);
  }

  /** 
   *
   */
  virtual double compute( const size_t varA, const size_t varB ) = 0;

  /** Return the number of total elements
   *
   */
  virtual size_t size() const = 0;

  /** Invalides any current caching
   *
   */
  virtual void invalidCache() = 0;

};

/////////////////////////////////////////////////////

/** Common Interface for Distance Matrix
 *  Usage: DissimilarityMatrix(a,b): return either the similarity or dissimilarity
 *  between 2 objects indexed by a and b
 */
struct DissimilarityMatrix: public CompMatrix {};

/////////////////////////////////////////////////////

/** Common Interface for Distance Matrix
 *  Usage: DissimilarityMatrix(a,b): return either the similarity or dissimilarity
 *  between 2 objects indexed by a and b
 */
struct SimilarityMatrix: public CompMatrix {};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_COMPARABLE_HPP

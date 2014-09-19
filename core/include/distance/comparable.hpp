/****************************************************************************************
 * File: comparable.hpp
 * Description: provides the common interface for CompareMatrix
 * @author: siolag161 (thanh.phan@outlook.com) /* affiliation */
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_COMPARABLE_HPP
#define SAMOGWAS_COMPARABLE_HPP
 
#include "stddef.h" // size_t
namespace samogwas
{

/** Common Interface for Distance Matrix and Similarity Matrix.
 * CS shoulb be indicated too in the headfile
 *  Usage: CompMatrix(a,b): returns either the similarity or dissimilarity
 *  between two objects indexed by a and b
 */
struct CompMatrix { // bad identifier
  
  /** Computes and returns either the similarity or distance between two objects  
   * indexed by varA and varB in a given matrix. This is a proxy method, that in turn calls the compute method
   * CS the term proxy is not common to every potential reader. Brief explanation (and motivation) expected
   */
  virtual double operator()( const size_t varA, const size_t varB ) {
    return compute(varA, varB);
  }

  /** 
   *
   */
  virtual double compute( const size_t varA, const size_t varB ) = 0;

  /** Returns the number of total elements
   * CS vague
   *
   */
  virtual size_t size() const = 0;

  /** Invalides any current caching
   * CS not informative
   */
  virtual void invalidCache() = 0;

};

/////////////////////////////////////////////////////

/** Common Interface for Distance Matrix
 *  Usage: DissimilarityMatrix(a,b): returns the dissimilarity
 *  between two objects indexed by a and b
 */
struct DissimilarityMatrix: public CompMatrix {};

/////////////////////////////////////////////////////

/** Common Interface for Distance Matrix
 *  Usage: SimilarityMatrix(a,b): returns the similarity
 *  between two objects indexed by a and b
 */
struct SimilarityMatrix: public CompMatrix {};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_COMPARABLE_HPP

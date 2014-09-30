/****************************************************************************************
 * File: mutual_information_similarity.hpp
 * Description: Similarity based on mutual information metric. 
 *              If a positive threshold is specified, the similarity values are discrete and belong to {0,1}.
 *              Otherwise, the actual continuous similarity values are kept.
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 12/06/2014
 ***************************************************************************************/
#ifndef SAMOGWAS_MUTUAL_INFORMATION_SIMILARITY_HPP
#define SAMOGWAS_MUTUAL_INFORMATION_SIMILARITY_HPP

#include "information_dissimilarity.hpp" // MutInfoDissimilarity
#include "comparable.hpp" // SimilarityMatrix

namespace samogwas
{

// @todo: change: pos-> positions, maxPos -> maxDist
template<class DataMatrix>
struct MutInfoSimilarity: public SimilarityMatrix {
  MutInfoSimilarity( std::vector< std::vector<int> >& dm, std::vector<int>& pos, unsigned maxPos, double thres ):
      m_miDiss( dm, pos, maxPos, thres) {}

  virtual double compute( const size_t varA, const size_t varB ) {
    return 1.0 - m_miDiss(varA,varB);
  }

  virtual size_t nbrVariables() const { return m_miDiss.nbrVariables(); }
  virtual void invalidate() { m_miDiss.invalidate(); }
  
 private:
  MutInfoDissimilarity<DataMatrix> m_miDiss;
};

} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_MUTUAL_INFORMATION_SIMILARITY_HPP

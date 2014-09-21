/****************************************************************************************
 * File: information_similarity.hpp
 * Description: Information similarity provides the similarity based on mutual information.
 * -----------  It can be obtained via Information Dissimilarity as Simi(a,b) = 1.0 - Diss(a,b)
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 30/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_INFORMATION_SIMILARITY_HPP
#define SAMOGWAS_INFORMATION_SIMILARITY_HPP

#include "information_dissimilarity.hpp"
#include "comparable.hpp"
namespace samogwas
{

template<class DataMatrix>
struct MutInfoSimilarity: public SimilarityMatrix {
  MutInfoSimilarity( std::vector< std::vector<int> >& dm, std::vector<int>& pos, unsigned maxPos, double thres ):
      m_miDiss( dm, pos, maxPos, thres) {}

  virtual double compute( const size_t varA, const size_t varB ) {
    return 1.0 - m_miDiss(varA,varB);
  }

  virtual size_t size() const { return m_miDiss.size(); }
  virtual void invalidCache() { m_miDiss.invalidCache(); }
  
 private:
  MutInfoDissimilarity<DataMatrix> m_miDiss;
};

} // namespace samogwasends here. samogwas

/****************************************************************************************/
#endif // SAMOGWAS_INFORMATION_SIMILARITY_HPP

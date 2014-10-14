/****************************************************************************************
 * File: euclidian.hpp
 * Description: 
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 13/10/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_EUCLIDIAN_HPP
#define SAMOGWAS_EUCLIDIAN_HPP

#include <math.h> // std::sqrt

namespace samogwas
{

template<class DataMatrix>
struct EuclidianDissimilarity: public DissimilarityMatrix {
  EuclidianDissimilarity( DataMatrix& dm );

  virtual double compute( const size_t varA, const size_t varB );
  
  virtual size_t nbrVariables() const {
    return dataMat.size();
  }
  
  virtual void invalidate() {}


 private:
  // reference to the actual data
  DataMatrix& dataMat;
};

} // namespace samogwasends here. samogwas
///////////////////////////////////////////////////////////////////////

namespace samogwas
{

template<class DM>
EuclidianDissimilarity<DM>::EuclidianDissimilarity( DM& dm): dataMat(dm) {}


template<class DM>
double EuclidianDissimilarity<DM>::compute( const size_t varA, const size_t varB ) {
  double sumOfSquares = 0.0;

  size_t nbrCols = dataMat.at(varA).size();
  for ( size_t c = 0; c < nbrCols; ++c ) {
    sumOfSquares += (dataMat[varA][c] - dataMat[varB][c])*(dataMat[varA][c] - dataMat[varB][c]);
  }

  return std::sqrt(sumOfSquares);
}

} // namespace samogwas ends here.

/****************************************************************************************/
#endif // SAMOGWAS_EUCLIDIAN_HPP

/****************************************************************************************
 * File: mutual_information.hpp
 * Description: This module provides methods to compute the mutual information of two variables.
 * MI(X,Y) = H(X) + H(Y) - H(X,Y)
 * 
 * @author: @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 29/12/2013

 ***************************************************************************************/
#ifndef SAMOGWAS_MUTUAL_INFORMATION_HPP
#define SAMOGWAS_MUTUAL_INFORMATION_HPP //@todo: change the namespace

#include <math.h> //log
#include <map>
#include <stdlib.h> // abs
#include <algorithm>    // std::min
#include <numeric> // std::accumulate
#include <boost/shared_ptr.hpp>

#include "utils/type_utils.hpp" // for utility::Int2Type
#include "entropy.hpp"

namespace samogwas
{

template<int EstimateMethodType>
struct MutualInformation
{

  MutualInformation() {}

  /**
   * Takes two iterators as input and computes the mutual information
   * according to the method selected for estimating the entropy.
   */
  template<typename VIterator>
  double operator()( VIterator xBegin, VIterator xEnd, VIterator yBegin,  bool has_missing = false ) {
    return compute(xBegin, xEnd, yBegin, utility::Int2Type<EstimateMethodType>(), has_missing);
    
  }

  template<typename VecType>
  double operator()(const VecType& xVec, const VecType& yVec, bool has_missing = false )  {    
    return compute(xVec.begin(), xVec.end(), yVec.begin(), utility::Int2Type<EstimateMethodType>(), has_missing );
  }

  /** MatrixT passed as parameter is a row-major Matrix in which each Row denotes a variable.
  */
  template<template<class> class MatrixT, class T>
  boost::shared_ptr<MatrixT<double> > operator()(const MatrixT<T>&matrix, bool has_missing = false );
  
 protected:  
  template<typename VIterator>
  double compute( VIterator xBegin, VIterator xEnd, VIterator yBegin, utility::Int2Type<EMP>, bool has_missing = false);
  
  /** MatrixT passed as parameter is a row-major Matrix in which each Row denotes a variable.
  */
  template<template<class> class MatrixT, class T>
  boost::shared_ptr<MatrixT<double> > compute(const MatrixT<T>& mat, utility::Int2Type<EMP>, bool has_missing = false );
  
  //template<typename XIterator, typename YIterator>
  //double compute(XIterator xBegin, XIterator xEnd, YIterator yBegin, utility::Int2Type<DIRICHLET>);
  
  //template<template<class> class MatrixT, class T>
  //boost::shared_ptr<MatrixT<double> > compute(const MatrixT<T>& mat, utility::Int2Type<SCALED_MI>);

};

} // namespace samogwas ends here.

/****************************** IMPLEMENTATION BELOW THIS POINT **************************/
namespace samogwas
{

template<int EstimateMethodType> 
template<template<class> class MatrixT, class T>
boost::shared_ptr<MatrixT<double> > MutualInformation<EstimateMethodType>::operator()(const MatrixT<T>& mat, bool has_missing)
{
  return compute(mat, utility::Int2Type<EstimateMethodType>(), has_missing);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
template<int EstimateMethodType> 
template<typename VIterator>
double MutualInformation<EstimateMethodType>::compute( VIterator xBegin, VIterator xEnd,
                                                       VIterator yBegin, utility::Int2Type<EMP>,
                                                       bool has_missing )
{ 

  std::map<int, unsigned> xCountMap;
  std::map<int, unsigned> yCountMap;
  typedef std::pair<int,int> IntPair;
  std::map<IntPair, unsigned> jointCountMap;
  
  double vecLen = 0.0;
  for (; xBegin != xEnd; ++xBegin)  {
    const int xVal = *xBegin, yVal = *yBegin;
    if (xVal != DATA_MISSING_VALUE && yVal != DATA_MISSING_VALUE) {
      updateCountMap(xCountMap, xVal);
      updateCountMap(yCountMap, yVal);
      const IntPair jointKey = std::make_pair(*xBegin, *yBegin); 
      updateCountMap(jointCountMap, jointKey);
      vecLen++; 
    }
    yBegin++;      
  }
    
  double xLogSum = accumulate( xCountMap.begin(), xCountMap.end(), 0.0, sumLogCount<unsigned>);
  double yLogSum = accumulate( yCountMap.begin(), yCountMap.end(), 0.0, sumLogCount<unsigned>);
  double jointLogSum = accumulate( jointCountMap.begin(), jointCountMap.end(),
                                   0.0, sumLogCount<IntPair>);

  // double entropyX = log(vecLen) - 1/vecLen*(xLogSum);
  // double entropyY = log(vecLen) - 1/vecLen*(yLogSum);
  // double jEntropy = log(vecLen) - 1/vecLen*(jointLogSum);
  return log(vecLen) + 1/vecLen*(jointLogSum - xLogSum - yLogSum);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
template<int EstimateMethodType> 
template<template<class> class MatrixT, class T>
boost::shared_ptr<MatrixT<double> >
MutualInformation<EstimateMethodType>::compute( const MatrixT<T>& mat,
                                                utility::Int2Type<EMP>,
                                                bool has_missing )
{
  boost::shared_ptr<MatrixT<double> > result(new MatrixT<double>(mat.nbrRows(), mat.nbrRows()));
  std::map<int, double> entropyMap;

  Entropy<EMP> entropy;
  JointEntropy<EMP> mutualEntropy;
  
  for (unsigned i = 0; i < mat.nbrRows(); ++i) {    
    entropyMap[i] = entropy(mat[i], has_missing); 
  }
  
  for (unsigned i = 0; i < mat.nbrRows(); ++i) {  
    (*result)[i][i] = entropyMap[i]; // MI(X,X) = E(X)
    for (unsigned j = i+1; j < mat.nbrRows(); ++j) {
      double ijEntropy = mutualEntropy(mat[i], mat[j], has_missing);
       double mutInfo = entropyMap[i] + entropyMap[j] - ijEntropy;
       (*result)[i][j] = mutInfo; 
       (*result)[j][i] = (*result)[i][j];
    }
  }
  return result;
}

} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // MUTUAL_INFORMATION_HPP

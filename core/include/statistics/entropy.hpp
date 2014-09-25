/****************************************************************************************
 * File: entropy.hpp
 * Description: This module provides methods to compute the entropy of a given variable, and the joint
 * -------------entropy of two variables.
 * H(X) = - sum_i p(X = i) log (p(X = i)) 
 * H(X,Y) = - sum_i,j p(X = i, Y = j) log(p(X = i, Y = j))
 * 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 30/12/2013
 ***************************************************************************************/
#ifndef SAMOGWAS_ENTROPY_HPP
#define SAMOGWAS_ENTROPY_HPP

/** A set of methods for computing estimated entropies for random variables.
 * Method type is decided at compile-time via the template parameter.
 * Formula for computing Entropy: -sum(p_i * log_2(p_i)) 
 */

#include <math.h> // log
#include <map>
#include <stdlib.h> // abs
#include <algorithm>  // std::min
#include <numeric> // std::accumulate

#include "utils/type_utils.hpp" // utility::Int2Type

namespace samogwas
{

enum EstimationMethod {EMP = 0, DIRICHLET, SCALED_MI}; //EMP = empirical

/** Computes entropy in base 2 (@todo: log2).
 *
 */
template<int EstimationMethodT>
struct Entropy
{
  Entropy() {}

  template<typename XIterator>
  double operator()(XIterator xBegin, XIterator xEnd) {
    return compute(xBegin, xEnd, utility::Int2Type<EstimationMethodT>());
  }

  template<typename VecT>
  double operator()(const VecT& xVec) {
    return compute(xVec.begin(), xVec.end(), utility::Int2Type<EstimationMethodT>());
  }

 protected:  
  template<typename XIterator>
  double compute(XIterator xBegin, XIterator xEnd, utility::Int2Type<EMP>);
  
  template<typename XIterator>
  double compute(XIterator xBegin, XIterator xEnd, utility::Int2Type<DIRICHLET>);
  

};


/**
 * Formula for computing joint entropy of X and Y: -sum( p_xy( log_2(p_xy) ) )
 */
template<int EstimationMethodT>
struct JointEntropy
{
  JointEntropy() {}

  // We suppose that X and Y have the same size.
  template<typename XIterator, typename YIterator>
  double operator()(XIterator xBegin, XIterator xEnd, YIterator yBegin)
    { return compute(xBegin, xEnd, yBegin, utility::Int2Type<EstimationMethodT>()); }

  template<typename VecXT, typename VecYT>
  double operator()(const VecXT& xVec, const VecYT& yVec)
    { return (*this)(xVec.begin(), xVec.end(), yVec.begin()); }

 protected:  
  template<typename XIterator, typename YIterator>
  double compute(XIterator xBegin, XIterator xEnd, YIterator yBegin, utility::Int2Type<EMP>);
};


template<typename T>
void updateCountMap(std::map<T, unsigned>& countMap, const T& val);

template<typename T>
double sumLogCount(double &total, const std::pair<T, unsigned>& data);

} // namespace samogwas ends here. 

/****************************** IMLEMENTATION BELOW THIS POINT **************************/
namespace samogwas
{


template<int EstimationMethodT>
template<typename XIterator>
double Entropy<EstimationMethodT>::compute(XIterator xBegin, XIterator xEnd, utility::Int2Type<EMP>)
{
  std::map<unsigned, unsigned> xCountMap; 
  double vecLen = 0.0;
  for (; xBegin != xEnd; ++xBegin)
  {
    const unsigned xVal = *xBegin; 
    updateCountMap(xCountMap, xVal); 
    vecLen++;
  }
    
  double xLogSum = accumulate( xCountMap.begin(), xCountMap.end(), 0.0, sumLogCount<unsigned> );
  double enp = log(vecLen) - 1/vecLen*xLogSum;

  return enp;
}
/////////////////////////////////////////////////////////////////////////////

template<int EstimationMethodT>
template<typename XIterator, typename YIterator>
double JointEntropy<EstimationMethodT>::compute( XIterator xBegin, XIterator xEnd,
                                                 YIterator yBegin, utility::Int2Type<EMP>)
{

  typedef std::pair<unsigned, unsigned> IntPair;
  std::map<IntPair, unsigned> jointCountMap;
  
  double vecLen = 0.0;
  for (; xBegin != xEnd; ++xBegin)  {

    const IntPair jointKey = std::make_pair(*xBegin, *yBegin); 
    updateCountMap(jointCountMap, jointKey);
    vecLen++; yBegin++;      
  }
    
  double jointLogSum = accumulate(jointCountMap.begin(), jointCountMap.end(),
                                  0.0, sumLogCount<IntPair>);
   
  return log(vecLen) - 1/vecLen*(jointLogSum);
}


/////////////////////////////////////////////////////////////////////////////


/** @todo: change to find and 
 *
 */
template<typename T>
void updateCountMap(std::map<T, unsigned>& countMap, const T& val)
{
  if(!countMap.insert(std::make_pair(val, 1)).second) // if already in the map
  {
    ++(countMap)[val]; // inscrease the counter by 1 (otherwise put it into the map, set the counter to 1)
  }
}

/** @param counter is the type of the pointee of the countMap iterator. This pointee is a
 * pair(element, number of occurences of element)
 */
template<typename T>
double sumLogCount(double &total, const std::pair<T, unsigned>& counter)
{
  double count = double(counter.second);
  return total + count*log(count);
}


} // namespace samogwas ends here. 

/****************************************************************************************/
#endif // SAMOGWAS_ENTROPY_HPP

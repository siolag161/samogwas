/****************************************************************************************
 * File: latent_var_criteria.hpp
 * Description: This module provides different criteria for evaluating a latent variable inferred and imputed by FLTM.
 *
 * @author: Phan Duc Thanh (duc-thanh.phan@univ-nantes.fr) - Under supervision of Christine Sinoquet
 * @date: 06/01/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_LATENT_VAR_CRITERIA_HPP
#define SAMOGWAS_LATENT_VAR_CRITERIA_HPP

#include "statistics/entropy.hpp"
#include "statistics/mutual_information.hpp"

namespace samogwas
{

/** This functor computes (estimates) the average mutual information between the latent variable
 * and its children. The functor simply returns sum(scaledMutInfo) / nbrOfChildren.
 *
 */
struct AverageMutInfo {
/** This functor requiries as paratemers:
 *  - the vector of values imputed for the latent variable,
 *  - the data matrix,
 *  - a vector containing the children's indices (w.r.t the dataset).
 */
  template<class VecType, class MatrixType, class ClusterType>
  double operator()( const VecType &latentVar,
                     const MatrixType &mat,
                     ClusterType &childrenIdx );
};

} // namespace samogwas ends here.

/****************************** IMPLEMENTATION BELOW THIS POINT **************************/
namespace samogwas
{

/**
 *
 */
template<class VecType, class MatrixType, class ClusterType>
double AverageMutInfo::operator()( const VecType &latentVar,
                                   const MatrixType &mat,
                                   ClusterType &childrenIdx)
{
  if (childrenIdx.size() == 0) return 0.0;

  Entropy<EMP> entropy;
  JointEntropy<EMP> jointEntropy;
  MutualInformation<EMP> mutualInfo;
  
  double totalScaledMutInfo = 0.0;
  double varEntropy = entropy(latentVar);

  for ( auto& idx: childrenIdx) {
    double idxEntropy = entropy( mat[idx] );     
    double minEntropy = std::min( idxEntropy, varEntropy );

    if (minEntropy == 0)
      continue;

    double jEntropy = jointEntropy( latentVar, mat[idx] );
    double mutInfo = varEntropy + idxEntropy - jEntropy;
    double scaledMutualInfo = mutInfo / minEntropy;
    totalScaledMutInfo += scaledMutualInfo;
  }

  double result = totalScaledMutInfo / childrenIdx.size();
  return result;
}


} // namespace samogwas ends here.

/****************************************************************************************/
#endif // SAMOGWAS_LATENT_VAR_CRITERIA_HPP

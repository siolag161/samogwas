/****************************************************************************************
 * File: Distance.hpp // CS Happy to see this for the information_dissimilarity.hpp file
 * Description: Dissimilarity based on mutual information metric.
 * @author: siolag161 (thanh.phan@outlook.com)
 * @date: 12/06/2014
 * // CS specify that this dissimilarity is thresholded (< median -> 0, > median -> 1) (ex CAST_bin)
 * // specify for = median
 * // Is there some code elsewhere for unthresholded dissimilarities? (ex CAST_real)
 *

 ***************************************************************************************/
#ifndef CLUSTERING_DISTANCE_HPP // CS Yet another name!
#define CLUSTERING_DISTANCE_HPP

#include <stdlib.h> /* abs */
#include <map>

#include <boost/accumulators/accumulators.hpp> /** to compute median **/
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp> /* to compute median */

#include "comparable.hpp"
#include "statistics/mutual_information.hpp"

namespace samogwas
{

const static double MAX_DISTANCE = 1.0, MIN_DISTANCE = 0.0;

/** Metric values are scaled to the [0,1] range. 
 *  This method also takes into consideration the position ( physical distance in the genomics context )
 *  for computing. Normal usage, without this contraint, can be done by setting the constraints correspondingly,
 *  for example by setting the threshold as infinity. 
 * CS We can override this class to provide more flexibility also. // CS Is it done?
 */
template<class DataMatrix>
struct MutInfoDissimilarity: public DissimilarityMatrix {
  
  /** The construction which takes a reference to the actual dataset, the position.
   *  The thres value determines whether this returns a 0/1 distance or the actual distance.
   *  The 0/1 distance is obtainable if `thres` is set as a positive value and return `1` if actual
   *  value is greater than this threshold, otherwise returns 0. If `thres` negative then return actual
   *  distance.
   */
  MutInfoDissimilarity( DataMatrix& dm, std::vector<int>& pos, unsigned maxPos, double thres );

  
  /** Computes and returns the distance between 2 variables, indexed by varA and varB
   *  in the dataset.
   */
  virtual double compute( const size_t varA, const size_t varB );
  
  /** Returns the number of variables
   *
   */
  virtual size_t size() const { return dataMat.size(); }

  /** Invalidaes current caching values
   *
   */
  virtual void invalidCache() {
    distCache.clear();
    entropyMap.clear();
    entropyMap.resize(positions.size(), -MAX_DISTANCE);
    distCache = std::map< size_t, double >();
  }
  
 public:

  // reference to the actual data
  DataMatrix& dataMat;

  // reference to the positions data
  std::vector<int>& positions;

  // the maximum threshold for position
  unsigned maxPosition;

  // threshold for 0/1 data casting
  double m_thres;

  // for caching the entropy during computation
  std::vector<double> entropyMap;

  // for caching the distance during computation
  std::map< size_t, double > distCache;

  // the current median value
  double m_median;
};

  
template<class DM>
double mutualInformationDistance( std::vector<double>& entropyMap,
                                  std::map<size_t,double>& distCache,
                                  const DM& dataMat,
                                  const size_t varA,
                                  const size_t varB );
template<class DM>
MutInfoDissimilarity<DM>::MutInfoDissimilarity( DM& dm,
                                                std::vector<int>& pos,
                                                unsigned maxPos,
                                                double thres ):
    dataMat(dm), positions(pos),
    maxPosition(maxPos), entropyMap(std::vector<double> (pos.size(), -MAX_DISTANCE)),
    m_thres(thres)
{
  size_t nbrVars = positions.size();
  if ( thres > 0 ) { // binary-->median, have to pre-compute all    
    boost::accumulators::accumulator_set< double,
                                          boost::accumulators::stats<boost::accumulators::tag::median(boost::accumulators::with_p_square_quantile) > > acc;
#pragma omp parallel for
    for ( size_t varA = 0; varA < nbrVars; ++varA) {
      for ( size_t varB = varA+1; varB < nbrVars; ++varB ) {
        if ( abs( positions[varA] - positions[varB] ) > maxPos ) continue; 
        double result = mutualInformationDistance( entropyMap, distCache, dataMat, varA, varB );
#pragma omp critical
        acc(result);
      }
    }
#pragma omp critical
    m_thres = boost::accumulators::median(acc);
  }
}

template<class DM>
double MutInfoDissimilarity<DM>::compute( const size_t varA, const size_t varB ) {
  if ( varA == varB ) return MIN_DISTANCE;
  if( varA > varB ) return this->compute( varB, varA );
  if ( abs( positions[varA] - positions[varB] ) >= maxPosition ) return MAX_DISTANCE;

  size_t nbrVars = dataMat.size();
  size_t key = 2*nbrVars*varA + varB;   
  double result = MAX_DISTANCE, rs = MAX_DISTANCE;
  std::map<size_t, double>::iterator iter;
      
  // #pragma omp critical
  iter = distCache.find(key);

  if ( iter == distCache.end() ) {
    rs =  mutualInformationDistance( entropyMap, distCache, dataMat, varA, varB );
  }

  // #pragma omp critical
  distCache[key] = rs;
    
  if ( m_thres > 0 ) { // binary-->median, have to pre-compute all
    result = ( rs > m_thres ) ? MAX_DISTANCE : MIN_DISTANCE;  
  } else { // continue
    result = rs;   
  }

  return result;
}

  
template<class DM>
double mutualInformationDistance( std::vector<double>& entropyMap,
                                  std::map<size_t,double>& distCache,
                                  const DM& dataMat,
                                  const size_t varA,
                                  const size_t varB )
{  
  double result = samogwas::MAX_DISTANCE;
  size_t nVars = entropyMap.size();
  size_t key = 2*nVars*varA + varB;
  Entropy<EMP> entropy; // empirical entropy (estimates the probability by frequency)
  JointEntropy<EMP> mutEntropy; // empirial mutual entropy;
  double enA, enB;
#pragma omp critical
  {  
    if ( entropyMap.at(varA) < 0.0 ) {
      entropyMap[varA] = entropy( dataMat.at(varA) ); // computes entropy of varA only if not already done
    }
    enA = entropyMap[varA];
  }

#pragma omp critical
  { 
    if ( entropyMap.at(varB) < 0.0 ) { // computes entropy of varB only if not already done
      entropyMap[varB] = entropy( dataMat.at(varB) ) ; // since this operation could be expensive
    }
    enB = entropyMap[varB];
  }

  double minEntropyAB = std::min(enA, enB); // takes the min
  if (minEntropyAB != 0) {
    double mutEntropyAB = mutEntropy( dataMat.at(varA), dataMat.at(varB) );
    double mutualInfoAB = enA + enB - mutEntropyAB; // classic formula
    double normalizedMutInfo = mutualInfoAB / minEntropyAB;        
    result = MAX_DISTANCE - normalizedMutInfo;
  }
  
  return result;
}



} // namespace clusteringends here. clustering

/****************************************************************************************/
#endif // CLUSTERING_DISTANCE_HPP

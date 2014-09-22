/****************************************************************************************
 * File: mutual_information_dissimilarity.hpp
 * Description: Dissimilarity based on mutual information metric. 
 *              If a positive threshold is specified, the dissimilarity values are discrete and belong to {0,1}.
 *              Otherwise, the actual continous dissimilarity values are kept.
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 12/06/2014
 ***************************************************************************************/
#ifndef SAMOGWAS_MUTUAL_INFORMATION_DISSIMILARITY_HPP 
#define SAMOGWAS_MUTUAL_INFORMATION_DISSIMILARITY_HPP

#include <stdlib.h> /* abs */
#include <map>

#include <boost/accumulators/accumulators.hpp> /** to compute median **/
#include <boost/accumulators/statistics/stats.hpp>  /** to compute median **/
#include <boost/accumulators/statistics/median.hpp> /** to compute median **/

#include "comparable.hpp"
#include "statistics/mutual_information.hpp"

namespace samogwas
{
const static double MAX_DISTANCE = 1.0, MIN_DISTANCE = 0.0; //@todo: DISTANCE -> DISSIMILARITY

/** Mutual information values are scaled to the [0,1] range, to obtain a metric.
 *  This method also takes into consideration the positions of the variables along a line (example: genome order).
 *  The dissimilarity between two variables is set to MAX_DISTANCE if these variables are 
 *  located too far apart (above maxDist).
 *  Usage without this maxDist contraint can be achieved by setting the value of maxDist to +infinity. 
 */
template<class DataMatrix>
struct MutInfoDissimilarity: public DissimilarityMatrix {
  
  /** The constructor takes a reference to the actual dataset, the positions, tha maxDist constraint and the threshold. 
   *  The threshold determines whether to compute the binary 0/1 dissimilarity values or the actual 
   *  continous dissimilarity values.
   *  The thresholded 0/1 dissimilarity is specified if `thres` is set as a positive value. 
   *  The dissimilarity value is `1` if the actual dissimilarity value is greater than this threshold; 
   *  otherwise, the dissimilarity value is 0. 
   *  If `thres` is negative, then the dissimilarity value is the actual dissimilarity value.
   */
  MutInfoDissimilarity( DataMatrix& dm, std::vector<int>& positions, unsigned maxDist, double thres );

  
  /** Computes and returns the dissimilarity between two variables, indexed by varA and varB
   *  in the dataset.
   */
  virtual double compute( const size_t varA, const size_t varB );
  
  /** Returns the number of variables // CS vague
   *
   */
  virtual size_t size() const { return dataMat.size(); }

  /** Invalidates current caching values // Cs which means?
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

  // reference to the positions data // CS Which means?
  std::vector<int>& positions;

  // the maximum threshold for position // CS Is it the maximal distance (in bp???) below which two objects are allowed
                                        // in the same cluster? Is it delta_dist?
  unsigned maxPosition;

  // threshold for 0/1 data casting
  double m_thres; // CS Could this be _thres, more simply?

  // for caching the entropy during computation // CS What do you mean by caching?
  std::vector<double> entropyMap; // CS Why map?

  // for caching the distance <CS dissimilarity???> during computation ?????
  std::map< size_t, double > distCache; // explain key and value

  // the current median value // CS of what? 
  double m_median; // CS Could this be _median, more simply?
};

  
template<class DM>
double mutualInformationDistance( std::vector<double>& entropyMap,
                                  std::map<size_t,double>& distCache,
                                  const DM& dataMat,
                                  const size_t varA,
                                  const size_t varB );
template<class DM>
MutInfoDissimilarity<DM>::MutInfoDissimilarity( DM& dm,
                                                std::vector<int>& pos, // bad identifier -> positions
                                                unsigned maxPos, // bad identifier
                                                double thres ):
    dataMat(dm), positions(pos),
    maxPosition(maxPos), entropyMap(std::vector<double> (pos.size(), -MAX_DISTANCE)), // CS unclear
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



} // namespace samogwas ends here

/****************************************************************************************/
#endif // SAMOGWAS_MUTUAL_INFORMATION_DISSIMILARITY_HPP 

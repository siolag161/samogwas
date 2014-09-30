/****************************************************************************************
 * File: g2.hpp
 * Description: G2 is a statistic used in an alternative to the chi-squared goodness-of-fit test.
 * ------------ The two test statistics usually have similar values and both have approximate
 *------------- chi-squared distributions when the model under test is a correct description
 *------------- of the data. 
 * @ref: http://en.wikipedia.org/wiki/G-test
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 22/06/2014

 ***************************************************************************************/
#ifndef STATS_G2_HPP
#define STATS_G2_HPP

#include <cmath>
#include <boost/math/distributions/chi_squared.hpp>
#include <numeric> // accumulate
#include <omp.h>
namespace stats
{

struct TestG2 {
 /** This function returns the statistic given by the formula: G2 = 2*sum( O_j*ln(O_j/E_j) )                  
 *   where the observed cell frequencies are denoted by O1, O2,..., Ok and the
 *   corresponding expected cell frequencies by E1, E2,..., Ek, respectively, and
 *   ln is the natural logarithm (Napierian). 
 */
  /** @param obs Observed frequencies (or counts) in each category
   *  @param exp Expected frequencies (or counts) in each category
   */
  template<class VectorType>
  double statistic( const VectorType& obs, const VectorType& exp )  {
    double statistic = 0.0;
    for (unsigned sz = 0; sz < obs.size(); ++sz) {
      if ( obs[sz]*exp[sz] != 0.0 ) {
        statistic += 2*obs[sz] * std::log( obs[sz] / exp[sz] );
      }
    }
    return statistic;
  }

  /////////////////////////////////////////////////////////////////////////
  /** @param obs Observed frequencies (or counts) in each category
   *  @param exp Expected frequencies (or counts) in each category
   */
  template<class VectorType>
  double gtest( const VectorType& obs, const VectorType& exp ) {
    double stat = statistic(obs, exp);
    const unsigned degreeFreedom = obs.size();
    return p_value(stat, degreeFreedom);
  }

  /////////////////////////////////////////////////////////////////////////
  /** This function returns the p-value of observing a value of a random variable
   *  which follows the chi-squared distribution.
   */
  double p_value( const double statistic, const unsigned int degreeFreedom )  {
    boost::math::chi_squared dist(degreeFreedom);
    double p_value = 1.0;
    try {
      p_value = 1 - boost::math::cdf(dist, statistic);
    }
    catch ( std::exception& e) {
      p_value = 1.0;
    }
    return p_value;
  }

 ////////////////////////////////////////////////////////////////////////////////
  /** This function returns the statistic g2 = 2*sum( O_j*ln(O_j/E_j) )                  
 *   where the observed cell frequencies are denoted by O1, O2,..., Ok and the
 *   corresponding expected cell frequencies by E1, E2,..., Ek, respectively, and
 *   ln is the natural logarithm (Napierian). 
   *  @param contigencyTab 
   *  @param useYates boolean to indicate whether the Yates correction is used.
   *         The effect of Yates' correction is to prevent overestimation of statistical
   *         significance for small data (at least one cell of the table has an expected count smaller than 5).
   */
  template<class ContingencyTabT>
  double gtest( const ContingencyTabT& contingencyTab, bool useYates = false ) {
    assert(contingencyTab.size() > 0);
    const unsigned nbrRows = contingencyTab.size();
    const unsigned nbrColumns = contingencyTab[0].size();

      ContingencyTabT tab = contingencyTab;
    if ( useYates ) { // Yates' correction is only valid when we have a 2x2 contingency table.
      assert( nbrRows == 2 ); assert( nbrColumns == 2); 
      if ( tab[0][0]*tab[1][1] > tab[0][1]*tab[1][0] ) {
        tab[0][0] -= 0.5; tab[1][1] -= 0.5;
        tab[0][1] += 0.5; tab[1][0] += 0.5;
      } else {
        tab[0][0] += 0.5; tab[1][1] += 0.5;
        tab[0][1] -= 0.5; tab[1][0] -= 0.5;
      }
    }
    double tableSum = 0.0;
    std::vector<double> rowSums(nbrRows, 0.0), columnSums(nbrColumns, 0.0);

    for (unsigned row = 0; row < nbrRows; ++row) {
      for (unsigned col = 0; col < nbrColumns; ++col) {
        double caseVal = tab[row][col];
          columnSums[col] += caseVal;
          rowSums[row] += caseVal;
          tableSum += caseVal;        
      }
    }

    double statistic = 0.0;
    for (unsigned row = 0; row < nbrRows; ++row) {
      for (unsigned col = 0; col < nbrColumns; ++col) {
        double observed = tab[row][col];
        double expected = rowSums[row]*columnSums[col] / tableSum;
        if ( observed*expected != 0.0 ) {
          statistic += 2 * observed * std::log( observed / expected );
        }
      }
    }    
    const unsigned degreeFreedom = (nbrRows-1)*(nbrColumns-1);
    return p_value(statistic, degreeFreedom);
  }  
};

} // namespace stats ends here. 

/****************************************************************************************/
#endif // STATS_G2_HPP


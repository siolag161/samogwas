/****************************************************************************************
 * File: chisq.hpp
 * Description: Templated chi-squared contingency table tests 
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 24/06/2014

 ***************************************************************************************/
#ifndef STATS_CHISQ_HPP
#define STATS_CHISQ_HPP

#include <cmath>

namespace stats
{

struct TestChiSquared {

  /** This function returns the statistic  chi2 = sum_{i=1}^{n} (O_i - E_i)^2 / E_i, 
   *  where n is the number of cells in the contingency table. O_i is the ith observed value in vector obs, and 
   *  E_i is the ith expected value in exp vector.
   *  @param obs Observed frequencies (or counts) in each category
   *  @param exp fExpected frequencies (or counts) in each category
   *  @param useYates boolean to indicate whether the Yates correction is used.
   *         The effect of Yates' correction is to prevent overestimation of statistical
   *         significance for small data (at least one cell of the table has an expected count smaller than 5).
   */
  template<class VecType>
  double statistic( const VecType& obs, const VecType& exp, bool useYates = false  ) const {
    double statistic = 0.0;
    for (unsigned sz = 0; sz < obs.size(); ++sz) {
      if ( obs[sz]*exp[sz] != 0.0 ) {
        double t = 0.0;
        if (useYates)
          t = abs(obs[sz] - exp[sz]) - 0.5;
        else
          t = obs[sz] - exp[sz];
        statistic += t*t / exp[sz];
      }
    }
    return statistic;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  /** This function returns the p-value by computing the statistic which, under the null hypothesis, follows
   *  the chi-squared distribution.
   *  @param obs Observed frequencies (or counts) in each category
   *  @param exp Expected frequencies (or counts) in each category
   */
  template<class VecType>
  double chisqTest( const VecType& obs, const VecType& exp, bool useYates = false ) const {
    double stat = statistic(obs, exp, useYates);
    const unsigned degreeFreedom = obs.size();
    return p_value(stat, degreeFreedom);
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  /** This function returns the p-value of observing a value of a random variable
   *  which follows the chi-squared distribution.
   */
  double p_value( const double statistic, const unsigned int degreeFreedom ) const {
    boost::math::chi_squared dist(degreeFreedom);
    double p_value = 1.0;
    try {
      p_value = 1-boost::math::cdf(dist, statistic);
    }
    catch ( std::exception& e) {
      p_value = 1.0;
    }
    return p_value;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  /** This function returns the statistic  chi2 = sum_{i=1}^{n} (O_i - E_i)^2 / E_i, 
   *  where n is the number of cells in the contingency table. O_i is the observed value and 
   *  E_i is the expected value.
   *  @param contigencyTab 
   *  @param useYates boolean to indicate whether the Yates correction is used.
   *         The effect of Yates' correction is to prevent overestimation of statistical
   *         significance for small data (at least one cell of the table has an expected count smaller than 5).
   */
  template<class ContingencyTabT>
  double chisqTest( const ContingencyTabT& contingencyTab, bool useYates = false ) const {
    assert(contigencyTab.size() > 0);
    const unsigned nbrRows = contingencyTab.size();
    const unsigned nbrColumns = contigencyTab[0].size();

    ContigencyTabT tab = contigencyTab;    
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
        if ( expected != 0.0 ) { 
          double t = 0.0;
          if (useYates)
            t = std::abs(observed - expected) - 0.5;          
          else
            t = observed - expected;
          statistic += t*t / expected;
        }        
      }
    }

    
    const unsigned degreeFreedom = (nbrRows-1)*(nbrColumns-1);
    return p_value( statistic, degreeFreedom );
  }  
};

} 

/****************************************************************************************/
#endif // STATS_CHISQ_HPP

/****************************************************************************************
 * File: test_statistics.hpp
 * Description: This module provides the common interface for all the available templated association tests:
 * -----------  Fisher exact test, ChiSquared (with and without Yates' correction), 
 * -----------  G2 (with and without Yates' correction).
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 22/06/2014

 ***************************************************************************************/
#ifndef STATS_TEST_STATISTICS_HPP
#define STATS_TEST_STATISTICS_HPP

#include <vector>
#include "g2.hpp"
#include "fisher.hpp"
#include "chisq.hpp"

namespace stats
{

template <int I>
struct Int2Type
{
  enum { value = I };
};

enum TestType { G2 = 0, G2_YATES, CHISQ, CHISQ_YATES, FISHER };

template< int Test >
struct StatisticTest {
  
  static const int test = (Test == G2) ? G2 :
                          (Test == G2_YATES) ? G2_YATES :
                          (Test == CHISQ) ? CHISQ :
                          (Test == CHISQ_YATES) ? CHISQ_YATES : FISHER;

  template<class ContingencyTableType>
  double operator()( ContingencyTableType& mat ) const;

  /** @param
   *  @param 
   * requires that data.size() == which.size()
   */ 
  template<class VecType>
  double operator()( const VecType& geno,
                     const VecType& pheno,
                     const int cardGenotype,
                     const int cardPhenotype ) const {
    std::vector< std::vector<double> > contingencyTab( cardGenotype, 
                                                       std::vector<double>(cardPhenotype, 0.0) );
    for (int i = 0; i < geno.size(); ++i) {
      int row = geno[i];
      int col = pheno[i];
      ++contingencyTab[row][col];
    };
    return p_value( contingencyTab, Int2Type<test>() );    
  }
  
 private:  
  template<class ContingencyTableType>
  double p_value( ContingencyTableType& mat, Int2Type<G2> ) const;

  template<class ContingencyTableType>
  double p_value( ContingencyTableType& mat, Int2Type<G2_YATES> ) const;
  
  template<class ContingencyTableType>
  double p_value( ContingencyTableType& mat, Int2Type<CHISQ> ) const;

  template<class ContingencyTableType>
  double p_value( ContingencyTableType& mat, Int2Type<CHISQ_YATES> ) const;
  
  template<class ContingencyTableType>
  double p_value( ContingencyTableType& mat, Int2Type<FISHER> ) const;

};

} // namespace stats ends here. 

/****************************** IMLEMENTATION BELOW THIS POINT **************************/
namespace stats
{

template< int Test >
template<class ContingencyTableType>
double StatisticTest<Test>::p_value( ContingencyTableType& mat, Int2Type<G2> ) const {
  TestG2 g2;
  return g2.gtest(mat, false);
}


template< int Test >
template<class ContingencyTableType>
double StatisticTest<Test>::p_value( ContingencyTableType& mat, Int2Type<G2_YATES> ) const {
  TestG2 g2;
  return g2.gtest(mat, true);
}

template< int Test >
template<class ContingencyTableType>
double StatisticTest<Test>::p_value( ContingencyTableType& mat, Int2Type<CHISQ> ) const {
  TestChiSquared chisq;
  return chisq.chisqTest(mat, false);
}


template< int Test >
template<class ContingencyTableType>
double StatisticTest<Test>::p_value( ContingencyTableType& mat, Int2Type<CHISQ_YATES> ) const {
  TestChiSquared chisq;
  return chisq.chisqTest(mat, true);
}

template< int Test >
template<class ContingencyTableType>
double StatisticTest<Test>::p_value( ContingencyTableType& mat, Int2Type<FISHER> ) const {
  TestFisher fisher;
  return fisher.fisherTest(mat);
}


} // namespace stats ends here. 

/****************************************************************************************/
#endif // STATS_TEST_STATISTICS_HPP


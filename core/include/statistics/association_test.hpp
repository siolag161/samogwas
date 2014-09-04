/****************************************************************************************
 * File: association_test.hpp
 * Description: This modules contains the following statistical of genotype association 
 * @author: Phan Duc Thanh (duc-thanh.phan@univ-nantes.fr) - Under supervision of Christine Sinoquet (christine.sinoquet@univ-nantes.fr)
 * @date: 25/06/2014

 ***************************************************************************************/
#ifndef STATS_ASSOCIATION_TEST_HPP
#define STATS_ASSOCIATION_TEST_HPP

#include "statistics.hpp"

namespace stats
{

struct StatTest {
  
  typedef std::vector<int> PhenoVec;
  typedef std::vector<int> GenoVec;
  typedef std::vector< GenoVec > GenoMat;

 public:
  
  /**
   *
   */
  StatTest() { }

  /**
   *
   */
  ~StatTest() {}

  /** Takes 2 vectors as input and performs a contingency table test. 
   *  Returns the p-value
   */
  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) const = 0;

  /**
   *
   */
  virtual double operator()( const GenoVec& geno,
                             const PhenoVec& pheno,
                             const unsigned cardGenotype,
                             const unsigned cardPhenotype ) const {
    return execute(geno, pheno, cardGenotype, cardPhenotype);
  }

 public:
  std::string name;
};


/** Performs a chi-squared contingency table test for count data
 *
 */
struct ChiSq: public StatTest  {
  
  ChiSq() { name = "ChiSq"; }
  
  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) const {
    return chisq(geno, pheno, cardGenotype, cardPhenotype);
  } 
  
  stats::StatisticTest<stats::CHISQ> chisq;
};


/** Performs a chi-squared contingency table test for count data, with Yates-correction
 *
 */
struct ChiSqCor: public StatTest  {
  
  ChiSqCor() { name = "ChiSq-Yates"; }  
  
  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) const {
    return chisq(geno, pheno, cardGenotype, cardPhenotype);
  }
  
  stats::StatisticTest<stats::CHISQ_YATES> chisq;
};

/** Performs a Fisher exact test for count data
 *  The underlying implementation borrows code from the R project
 */
struct Fisher: public StatTest  {
  Fisher() { name = "Fisher"; }

  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) const { return fisher(geno, pheno, cardGenotype, cardPhenotype); }
 
  stats::StatisticTest<stats::FISHER> fisher;
};

/** Performs a G-test for count data
 */
struct G: public StatTest  {
  G() { name = "G2"; }
  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) const { return g2(geno, pheno, cardGenotype, cardPhenotype); }
  
  stats::StatisticTest<stats::G2> g2;
};

/** Performs a G-test for count data, with Yates correction
 */
struct G2Cor: public StatTest  {
  
  G2Cor() { name = "G2-Yates"; }
  
  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) {
    return g2( geno, pheno, cardGenotype, cardPhenotype );
  }
  stats::StatisticTest<stats::G2_YATES> g2;
};

} // namespace stats ends here. stats


/****************************************************************************************/
#endif // STATS_ASSOCIATION_TEST_HPP

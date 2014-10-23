/****************************************************************************************
 * File: association_test.hpp
 * Description: This module contains the different tests for genotype-phenotype associations/ 
 * @author: Adapted by Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 25/06/2014

 ***************************************************************************************/
#ifndef STATS_ASSOCIATION_TEST_HPP
#define STATS_ASSOCIATION_TEST_HPP

#include "test_statistics.hpp"

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

  /** Takes 2 vectors as input and performs an independance test. 
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

///////////////////////////////////////////////////////////////////////////////////////
/** Performs a chi-squared test.
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

///////////////////////////////////////////////////////////////////////////////////////
/** Performs a chi-squared with Yates' correction.
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

///////////////////////////////////////////////////////////////////////////////////////
/** Performs a Fisher exact test.
 */
struct Fisher: public StatTest  {
  Fisher() { name = "Fisher"; }

  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) const {
    return fisher(geno, pheno, cardGenotype, cardPhenotype);
  }
 
  stats::StatisticTest<stats::FISHER> fisher;
};

///////////////////////////////////////////////////////////////////////////////////////
/** Performs a G-test.
 */
struct G: public StatTest  {
  G() { name = "G2"; }
  virtual double execute( const GenoVec& geno,
                          const PhenoVec& pheno,
                          const unsigned cardGenotype,
                          const unsigned cardPhenotype ) const { 
    return g2(geno, pheno, cardGenotype, cardPhenotype);
  }
  
  stats::StatisticTest<stats::G2> g2;
};

///////////////////////////////////////////////////////////////////////////////////////
/** Performs a G-test with Yate' correction.
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

////////////////////////////////////////////////////////////////////

struct GWASAssociationTest {

  typedef std::vector<int> PhenoVec;
  typedef std::vector<int> GenoVec;
  typedef std::vector< GenoVec > GenoMat;
  typedef std::vector<GenoVec> DataMat;

 public:
  GWASAssociationTest( const DataMat& data,
                       const PhenoVec& pheno,
                       const unsigned cardPheno ): dataMat(data),
                                                   phenoVec(pheno),
                                                   phenoCard(cardPheno) {}

  
  /** Takes 2 vectors as input and performs an independance test. 
   *  Returns the p-value
   */
  virtual double execute( const size_t& genoIdx, const unsigned genoCard ) const {
    
    return statTest->execute( dataMat[genoIdx], phenoVec, genoCard, phenoCard);
  }

 private:
  DataMat dataMat;
  std::shared_ptr<StatTest> statTest;
  const PhenoVec& phenoVec;
  unsigned phenoCard;
};


} // namespace stats ends here.


/****************************************************************************************/
#endif // STATS_ASSOCIATION_TEST_HPP

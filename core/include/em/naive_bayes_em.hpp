/****************************************************************************************
 * File: naive_bayes_em.hpp
 * Description: Specific instance of the EM algorithm dedicated to the parameter learning of
 * -----------the Naive Bayes multinomial model.
 *
 * @author: Duc-Thanh Phan siolag161 (thanh.phan@outlook.com), under the supervision of Christine Sinoquet
 * @date: 31/07/2014

 ***************************************************************************************/
#ifndef SAMOGWAS_NAIVE_BAYES_EM_HPP
#define SAMOGWAS_NAIVE_BAYES_EM_HPP

#include "core_em.hpp"
namespace samogwas
{


struct NaiveBayesEM: public EMInterface {
  typedef std::vector< std::vector<int> > Matrix;

 public:
  NaiveBayesEM( int nRestarts, int imMode ): nbrRestarts(nRestarts), imputMode(imMode) {}
  ~NaiveBayesEM() {}

 protected:
  virtual void run( ResultEM& result,
                    const Variable& latentVar,
                    const Variables& variables,
                    const Matrix& dataTable,
                    const double threshold,
                    std::vector< std::vector<bool> > *&defTable);


  virtual void imputeLatent( ResultEM& result,
                             const plSymbol& latentVar,
                             const Matrix& dataTable,
                             EMLearner& bestModel,
                             plMatrixDataDescriptor<int> &dataDesc );

  int nbrRestarts;
  int imputMode; // methods for imputing missing values  (ARGMAX or DRAW)

};

} // namespace samogwas ends here. 


/****************************************************************************************/
#endif // SAMOGWAS_NAIVE_BAYES_EM_HPP
